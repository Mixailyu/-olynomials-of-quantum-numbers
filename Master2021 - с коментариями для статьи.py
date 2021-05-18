#!/usr/bin/python
# -*- coding: utf-8 -*- # необходима для ввода коментариев на русском


## Сначала заведем класс Vector(n) - это объект содержащий в себе иформацию
## о неком невозмущенном векторе n, а именно:
## информацию о числах заполнения self._vec - характеризующих энергетическое состояние (n__1,n__2,...);
## информацию о постоянных множителях self._const - ангармонические постоянные, производные функции дипольного момента;
## информацию о числовых множителях self._NF возникших в результате воздействия на вектор операторов рождения и уничтожения;
## В общем это объект обладающий свойствами бра или кет вектора описывающего состояние системы в энергетическом представлении.
## В этот объект встроены функции которые могут вызвать информацию хранящююся внутри этого объекта: vec(self), NF(self), const(self)
class Vector():
  ''' Класс расчета NF -numeric factor, числовой множитель, после применения оператора к вектору '''
  def __init__(self, vec):
    self._vec = list(vec)
    self._NF=[]
    self._const=[]
  def vec(self):
    return tuple(self._vec)
  def NF(self):
    return '*'.join(['sqrt(n__%s)' % (VECTOR_INDEX[i] + num2str(j)) for i,j in self._NF])
  def const(self):
    if len(self._const)>0:
      return '*'.join(['(%s)' % (i) for i in self._const if len(i)>0])
    else:return ''
## Теперь необходимо создать функции описывающие работу операторов рождения born(vec,i) и  уничтожения dead(vec,i)
## в качестве аргументов функция принимает объект класса Vector(), в данном случае обозначенный как vec, и индекс числа в векторе состояния
## на которое нужно подействовать
## врезультате этого воздействия объект vec притерпевает изменение, а именно:
## число заполнения n__i увеличивается/уменьшается на единицу при воздействии оператора рождения/уничтожения
## одновременно с этим записывается числовой множитель n__i+1 при воздействии оператора рождения или n__i при воздействии оператора уничтожения в список self._NF   
def born(vec,i):
  ''' Оператор рождения '''
  while len(vec._vec)<i+1:
      vec._vec.append(0)
  vec._vec[i] +=1
  vec._NF.append((i,vec._vec[i]))
def dead(vec,i):
  ''' Оператор уничтожения '''
  while len(vec._vec)<i+1:
      vec._vec.append(0)
  vec._NF.append((i,vec._vec[i]))
  vec._vec[i] -=1
## Так возмущения в теории полиномов квантовых чисел представлены в виде разложения в ряд по степеням нормальных координат q,
## нам необходимо завести функцию ksi^n где ksi=(dead,born) - более удобная вибрационная переменная равная q*sqrt(2),
## n - степень в которую нужно возвести ksi
def ksi_polinom(n:int)->list:
  ''' Разложение ksi()^j в полином  '''
  ksi = [dead,born]
  return [x[::-1] for x in product(ksi,repeat=n)]  # декартово произведение [[z,y,x] for x in a for y in a for z in a ]  

## Как видно в теории полиномов квантовых чисел очень часто нужно решать задачу о вычислении матричного элемента <bra|fksi(ijk)|ket>
## где fksi - оператор ksi^n, ijk-список индексов на которые нужно подействовать
## Для решения этой задачи была разработана функция CALC(bra,fksi,ket,ijk) её на вход подаются:
##  бра и кет вектора,которые могут иметь вид класса - Vector, последовательности не изменяемых элементов - tuple  или в виде списка - list;
##  список операторов fksi=ksi_polinom(n) которые должны подействовать на вектор кет;
##  список индексов ijk вектора кет на которые должны действовать операторы fksi;
## Важное замечание длинна списка ijk должна равняться степени n.
def CALC(bra,fksi,ket,ijk):
    ''' Считает матричный элемент <bra|fksi(ijk)|ket> fksi - оператор ksi()^j, ijk-список индексов на которые нужно подействовать '''
    const=''
    NF=''
    if isinstance(ket, Vector):
        const+=ket.const()
        NF+=ket.NF()
        ket=ket._vec
    if isinstance(bra, Vector):
        const+='*'+bra.const()
        NF+='*'+bra.NF()
        bra=bra._vec
    if isinstance(bra, tuple):bra=list(bra)
    if isinstance(ket, tuple):ket=list(ket)
    result = [Vector(ket) for i in range(len(fksi))]
    for vec,operations in zip(result,fksi):
        for op,i in zip(operations,ijk):
            op(vec,i) # применить оператор (born/dead) к элементу номер i вектора vec
    RESULT=[]
    for i in result:
        if bra==i._vec: RESULT.append('*'.join([c for c in ['1',i.NF(),const,i.const(),NF,] if c!='']))
    return '(%s)'% ('+'.join([i for i in RESULT if i!=''])) if len(RESULT)>0 else 0

##################################################
   
## Теперь когда мы задали базовые функции описывающие вектор состояния, операторы рождения и уничтожения,
##  а также операцию вычисления матричного элемента состоящего из невозмущенных векторов мы можем перейти п программираванию
##  основных рекурентных уравнений теории полиномов квантовых чисел (*),(**).
##  А именно уравнениий (*)для расчета проправки к энергии E(n,alfa) и (**) поправки к вектору состояния системы |n,alfa>,
##  где n - вектор состояния (n__1,n__2,...),  alfa - степень возмущения.

##   В выражениях (*),(**) для расчета проправки к энергии E(n,alfa) и поправки к вектору состояния системы |n,alfa> введен особое правило суммирования
##   (p,betta,gamma)alfa для поправки к энергии и (p,q,betta,gamma,nu)alfa для поправки к вектору состояния, которое говорит что сумма элементов внутри
##   скобок должна равняться значению alfa. Поэтому для большего удобства работы с индексами p,q,betta,.. и т.д. мы заведем именнованные кортежи
##   Index и Index2
from collections import namedtuple
Index = namedtuple('Index','p q betta gamma nu')
Index2 = namedtuple('Index2','p betta gamma')

## Для расчета возмущенного вектора состояния |n,alfa> (**) заведем функцию Perturbation_Vector для краткости запишим как AV(n,alfa)
## Первый цикл данной функции запускает перебор по всем возможным комбинациям индексов (p,q,betta,gamma,nu)
## с учетом ограничений в виде того что число p обязательно должно быть больше нуля и число q - должно быть четным.
## после идет повторный вызов функции AV(n,gamma) что эквивалентно вектору |n,gamma> из уравнения (**),
## определяется список операторов fksi действующих на вектор |n,gamma> и постоянный множитель p/alfa
## Второй цикл запускает перебор по всем возможным значениям вектора m таким, что m неравно n.
## На этом этапе применяется правило отбора, расчитывается значение функции DELTA(m,n,q),
## а также определяется значение вектора <m,betta| с помощью вызова функции AV(m,betta)
## Третий цикл запускает перебор индексов на которые должны воздействовать операторы fksi в векторе |n,gamma>.
## в этом цикле происходит расчет матричного элемента <m,betta|fksi(ijk)|n,gamma>
## с помощью вызова функции calc(BRA,KET,fksi,IJK,'s') данная функция является аналогом функции CALC с одним единственным отличием,
## она может работать возмущенными векторами 
                    #расчет поправки к вектору состояния#
#calculation Amendment to the Vector (AV)
def AV(ket=[0,0,0,0],ind=0):
  '''Расчитывает вектор |ket,ind>'''
  if ind==0:return [[tuple(ket),'']]#Vector(ket)
  N=[]
  for i in [Index(*i) for i in product(range(0,ind+1), repeat=5) if sum(i)==ind and i[1]%2==0 and i[0]>0]:#(p,q,бетта,гамма,ню)
    KET=AV(ket,i.gamma)
    fksi=ksi_polinom(i.p+2)
    NF1='(%s/%s)'%(i.p,ind)
    for bra in element_index(3*(i.betta+i.gamma)+i.p+2,ket): # i[0]+2=степени в которую будем возводить ksi
      if list(bra)==list(ket):continue
      if pravilo_otbora(bra,ket,i.p+2,i.betta,i.gamma)==0:continue
      DELTA=delta(bra,ket,i.q)
      if DELTA==0:continue
      BRA=AV(bra,i.betta)
  ##  if i.betta!=i.nu: P=AV(bra,i.nu)
##    else:P=[k for k in BRA]
      for ijk in INDEX(i.p,ket):
        A='(A__%s)'%(ijk) # это ангармоническая постоянная
        IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
        NF2=calc(BRA,KET,fksi,IJK,'s')
        if NF2==0:continue
        const='*'.join([k for k in [NF1,A,DELTA,NF2]])
        PP=[[V[0],'*'.join([k for k in [const,V[1]]if k!=''])] for V in AV(bra,i.nu)]
        N+=PP
  return N

def calc(BRA,KET,fksi,IJK,flag='s'):
  '''Расчитывает матричный элемент вида <BRA,alfa|fksi(IJK)|KET,betta>
     flag='s' вернуть рассчитанные матричные элементы в виде строки 'M1+M2+...'
     flag='l' вернуть рассчитанные матричные элементы в виде списка [M1,M2,...]'''
  result=[]
  for i in BRA:
    for j in KET:
      A=CALC(i[0],fksi,j[0],IJK)
      if A==0:continue
      I='*'.join([k for k in [A,i[1],j[1]] if k!=''])
      if len(I)>0:result.append(I)
  if flag=='s': return '(%s)'%('+'.join([k for k in result if k!=''])) if len(result)>0 else 0
  if flag=='l': return [k for k in result if k!='']
##################################################
        #Расчет поправки к энергии
  # calculation Amendment to the Energy (AE)

def AE(vec=[0,0,0,0],ind=0,flag='s'):
  '''Расчитывает вектор |ket,ind>'''
  if ind%2==1:return 0
  N=[]
  X=[]
  for i in [Index2(*i) for i in product(range(0,ind+1), repeat=3) if sum(i)==ind and i[0]>0]:#(p,бетта,гамма)
    KET=AV(vec,i.gamma)
    BRA=AV(vec,i.betta)
    fksi=ksi_polinom(i.p+2)
    NF1='(%s/%s)'%(i.p,ind)
    E_ijk=[]
    for ijk in INDEX(i.p,vec):
      A='(A__%s)'%(ijk)
      IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
      E=calc(BRA,KET,fksi,IJK,'s')
      if E==0:continue
      E_ijk.append(str(sy.simplify(E+'*'+A)))
    E_ijk='+'.join(E_ijk)
    N.append('%s*%s'%(NF1,E_ijk))
  if flag=='s': return '+'.join(N)
  if flag=='l': return N
##################################################         
def pravilo_otbora(bra,ket,nksi,alfa,betta):
    k = sum(bra)-sum(ket)
    return 0 if (alfa+betta+nksi) % 2 != k % 2 else 1

def delta(bra,ket,q):
    if q==0:
        omega = ['(%s*omega__%s)' % (r-l,VECTOR_INDEX[i])  for i,(l,r) in enumerate(zip(bra,ket)) if l!=r]
        return '(1/(%s))' % '+'.join(omega) if omega else 0
    if q==2:
        A=delta(bra,ket,0)
        B='-'.join([AE(ket,2,'s'),AE(bra,2,'s')])
        return '*'.join([A,A,B]) if A else 0
    else: return 0
##################################################
## добавленные функции и библиотеки
from itertools import product, combinations_with_replacement
import numpy as np

def num2str(n):
  ''' Число в строку со знаком, ноль - пусто'''
  return '%+d' % n if n else ''

VECTOR_INDEX = 'ijklmnopq'
VECTOR_INDEX_MAP = {k:i for i,k in enumerate(VECTOR_INDEX)}
operator_index = lambda i: [ VECTOR_INDEX_MAP[x] for x in i]

def INDEX(n:int,vec)->list:
  ''' Индексы вектора vec, степени возмущения n, vec дб tuple для хэширования '''
  return [''.join(sorted(i)) for i in combinations_with_replacement(VECTOR_INDEX[:len(vec)], n+2) ]    

def element_index(n:int,vec):
  '''генерирует список с возможными значениями k для вектора m_i=n_i+k_i, k_i принимает значение [-n:n] и сумма модулей k_i   меньше или равна n
  пример: element_index(1,(0,0))== [(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)]'''
  ''' Индексы элементов к которым будут применяться операторы '''
  rng = range(-n,n+1)
  el = [i for i in product(rng, repeat=len(vec)) if sum(map(abs,i)) <=n]
  el = np.array(el) + vec
  return list(map(tuple,el))
###############################################
import sympy as sy
A__iiii, A__iiij, A__iiik, A__iiil, A__iijj, A__iijk, A__iijl, A__iikk, A__iikl, A__iill, A__ijjj, A__ijjk, A__ijjl, A__ijkk, A__ijkl, A__ijll, A__ikkk, A__ikkl, A__ikll, A__illl, A__jjjj, A__jjjk, A__jjjl, A__jjkk, A__jjkl, A__jjll, A__jkkk, A__jkkl, A__jkll, A__jlll, A__kkkk, A__kkkl, A__kkll, A__klll, A__llll = sy.symbols('A__iiii, A__iiij, A__iiik, A__iiil, A__iijj, A__iijk, A__iijl, A__iikk, A__iikl, A__iill, A__ijjj, A__ijjk, A__ijjl, A__ijkk, A__ijkl, A__ijll, A__ikkk, A__ikkl, A__ikll, A__illl, A__jjjj, A__jjjk, A__jjjl, A__jjkk, A__jjkl, A__jjll, A__jkkk, A__jkkl, A__jkll, A__jlll, A__kkkk, A__kkkl, A__kkll, A__klll, A__llll')
A__iii, A__iij, A__iik, A__iil, A__ijj, A__ijk, A__ijl, A__ikk, A__ikl, A__ill, A__jjj, A__jjk, A__jjl, A__jkk, A__jkl, A__jll, A__kkk, A__kkl, A__kll, A__lll = sy.symbols('A__iii, A__iij, A__iik, A__iil, A__ijj, A__ijk, A__ijl, A__ikk, A__ikl, A__ill, A__jjj, A__jjk, A__jjl, A__jkk, A__jkl, A__jll, A__kkk, A__kkl, A__kll, A__lll')
n__i, n__j, n__k, n__l=sy.symbols('n__i, n__j, n__k, n__l')
omega__i, omega__j, omega__k, omega__l=sy.symbols('omega__i, omega__j, omega__k, omega__l')
###############################################
# расчет матричного элемента функции дипольного момента
# calculation the Matrix Element of the Dipole Moment Function (MEDMF)
def MEDMF (bra,ket,p, flag='s'):
	L=[]
	for i in range(p+1):
		for j in range(p+1):
			if i+j<=p:
				for k in range(1,p+2):
					if (i+j)%2==0 and k%2!=0:L.append((i,k,j))
					if (i+j)%2!=0 and k%2==0:L.append((i,k,j))
	BRA=[AV(bra,i) for i in range(p+1)]
	KET=[AV(ket,i) for i in range(p+1)]
	print('End BRA,KET')
	Dipol=[]
	for i in L:
		print(i)
		D_ijk=[]
		fksi=ksi_polinom(i[1])
		for ijk in INDEX(i[1]-2,ket):
			A='(D__%s)'%(ijk)
			IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
			D=calc(BRA[i[0]],KET[i[2]],fksi,IJK,'s')
			if D==0:continue
			D_ijk.append(str(sy.simplify(D+'*'+A)))
		Dipol.append('+'.join(D_ijk))
	if flag=='s': return '+'.join(Dipol)
	if flag=='l': return Dipol
############################ TEST ##########################
##
##N02=AV((0,0,0,0),2) #|n;2>
##N12=AV((1,0,0,0),2) #|n+1;2>
##N11=AV((1,0,0,0),1) #|n+1;1>
##N01=AV((0,0,0,0),1) #|n;1>
##N00=AV((0,0,0,0),0) #|n;0>
##N10=AV((1,0,0,0),0) #|n+1;0>
##
##fksi=ksi_polinom(1)
##IJK=[0]
##M10_ksi_02=calc(N10,N02,fksi,IJK,'l') #<n+1;0|ksi(i)|n;2>
##M01_ksi_11=calc(N01,N11,fksi,IJK,'l') #<n;1|ksi(i)|n+1;1>
##M00_ksi_12=calc(N00,N12,fksi,IJK,'l') #<n;0|ksi(i)|n+1;2>
##M=M10_ksi_02+M01_ksi_11+M00_ksi_12
##MM=[str(sy.simplify(i)) for i in M]
##MMM=[sy.simplify('+'.join(i)) for i in MM] #<n+1;0|ksi(i)|n;2>+<n;1|ksi(i)|n+1;1>+<n;0|ksi(i)|n+1;2>
## Теперь с помощью функции .subs(), библиотека sympy, можно заменить все буквенные константы на числа, тем самым получить ответ.

###########################  План на будущее  ##################################
##для нахождения поправок к вектору состояния более высокого порядка (3,4, и т.д.) необходимо
##определить функцию E(ind,vec)- поправка к эрнергии - это сделаю в течении пары недель,
##может быть раньше так как структура этой функции будет схожа с PP2
##определить функцию delta(bra,ket,q) для произвольного q - тут надо подумать.
## zamena=[(n__i,0),(n__j,0),(n__k,0),(n__l,0),(omega__i, 991.95),(omega__j, 496.67),(omega__k, 962.33),(A__iii, -229.5),(A__iij, -48.3),(A__ijj, -18.4),(A__ikk, -276.4),(A__jjj, -85.1), (A__jkk, -55.2),(A__iiii, 35.1),(A__iiij, 22.8),(A__iijj, -3.2),(A__ijjj, 7.3),(A__jjjj, 19.6),(A__iikk, 68.4),(A__ijkk, 16.4),(A__jjkk, -13.6),(A__kkkk, 43.9)]
