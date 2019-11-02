# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 08:23:44 2019

@author: jackl
GOAL: 平时学习密码学过程中涉及到的一些功能在 python 上进行实现
积累 
"""
import math
import random
import hashlib # 哈希函数库
class discreteMathematic_0110:
    def __init__(self):
        print('Discrete Mathematic By 0110')
        return
    def gcd(self,a,b):
        ''' 最大公约数 greatest common divider
        
        求解a,b 最大公约数 '''
        if b == 0:
            return a
        
        return self.gcd(b,a%b)
    def lcm(self,a,b):
        ''' 最小公倍数 least common multipler
        
        求解整数 a,b 的最小公倍数'''
        return int(a*b/self.gcd(a,b))
        
    def exGCD(self,a,b):
        '''拓展欧几里得算法
           
        输入两个正整数a,b，求解 ax+by=gcd(a,b)=g; 
        返回 x,y,g; x,y为方程的解,g 为a,b的最大公约数
        
        算法复杂度为 O(log(b)) '''
        if b == 0:
            return 1,0,a
        
        x,y,g = self.exGCD(b,a%b)
        return y,x-a//b*y,g
    
    def getMulInverse(self,a,b):
        ''' 求解乘法逆元
        
        a,b互素，ak%b=1 ,求 a 的逆元 k
        
        本质上就是构建 ax+by=1 这个方程，然后借助拓展欧几里得算法进行求解
        显然 要求 ax%b=1 中 x，就等价于求解 (ax+by)%b=1 中的 x，就等价于求 ax+by=1中的 x
        '''
        if self.gcd(a,b)!=1:
            print('a,b不互素。模数为b的情况下,a 不存在逆元')
            return -1
        
        x,y,g = self.exGCD(a,b)
        return x
    
    def printSTH(self):# 寻求 k ，并验证必然存在 k，使得 i^b % k = i % k 恒成立
        k = 10
        b = 5
        for i in range(2,100):
            print(int(math.pow(i,b))%k,i%k)
            
    def CRT(self,m,r):
        '''中国剩余定理(Chinese Remainder Theorem)
        r: 余数值元组，m: 模数元组
        返回值：模数元组和余数元组对应的元素 n
        
        满足 n%m[i]=r[i] 恒成立
        
        目标: 求解满足上述关系的最小非负整数 n
        
        天干地支纪年法与该定理有共性
        在孙子方法上进一步借助拓展欧几里德求得逆元，更为高效
        '''
        M = 1
        for item in m:
            M *= item

        goal = 0
        for i in range(0,len(m)):
            x,y,g = self.exGCD(int(M/m[i]),m[i])
            goal += int(M/m[i])*x*r[i]
        # 通过取模的方法，确定最小的非负整数目标值
        return goal%M
    
    def quickReductionExp(self,m,e,n):
        '''快速二分求幂运算
        
        求解 a^b%c
        
        这是一个递归算法
        算法时间复杂度为 O(logb)'''
        if e == 0: return 1
        r = self.quickReductionExp(m,e//2,n)
        r = r*r%n
        if e&1: # 是奇数
            return r*m%n
        
        return r
    def dlogAttack(self,m,c,n):
        '''RSA中离散对数
        
        m^e%n = c
        求 e
        穷举攻击
        时间复杂度 O(n)
        '''
        r = m
        k = 1
        while 1:
            if r == c: break
            r = r*m%n
            k += 1
        return k
        
    def eccAdd(self,P,Q,a,b,p):
        ''' 椭圆曲线加法 P+Q
        
        a,b 为椭圆曲线参数,为整数
        P,Q 为二维坐标点
        p 为模数，是一个素数
        返回 R = P+Q
        计算 P+Q
        注意： 返回 (0,0) 表示计算结果为 幺元
        
        '''
        xP,yP = P[0],P[1]
        xQ,yQ = Q[0],Q[1]
        
        # 处理 “幺元”情况
        if P==(0,0): return Q
        if Q == (0,0): return P
        if xP == xQ and yP != yQ:
            return (0,0)
        
        # 处理 “非幺元”情况
        if P == Q:
            if 3*xP*xP+a < 0: h = p -(-(3*xP*xP+a))%p
            else: h = (3*xP*xP+a)%p
            g = 2*yP%p
            g = self.getMulInverse(g,p)
        else:
            if yP - yQ < 0: h = p-(yQ-yP)%p
            else: h = (yP-yQ)%p
            if xP - xQ < 0: g = p-(xQ-xP)%p
            else: g = (xP-xQ)%p
            g = self.getMulInverse(g,p)
        lam = h*g%p # lambda 斜率

        xR = (lam*lam-xP-xQ)%p
        yR = (lam*(xP-xR)-yP)%p
        
        return (xR,yR)
        
    def printEccAllPoints(self,P,a,b,p):
        ''' 不停计算 P 点自加, 直到出现循环
            a,b 为椭圆曲线参数 ， p 为模数
        '''
        Q = P
        i = 1
        while 1:
            if Q == (0,0):
                print(i,'幺元')
            else:print(i,Q)
            Q = self.eccAdd(P,Q,a,b,p)
            
            if Q == P: break
            
            i= i+1
            
    def eccMultiply(self,k,P,a,b,p):
        ''' 椭圆曲线乘法（快速降幂思想）
        计算 kP 并返回
        a,b 为椭圆参数
        p 为 模数
        借助递归实现椭圆曲线乘法中的“快速降幂”
        
        时间复杂度 O(logk)
        
        '''
        if (k == 1): return P
        if (k == 0): return (0,0)
        R = self.eccMultiply(k//2,P,a,b,p)
        R = self.eccAdd(R,R,a,b,p)
        if (k&1): R = self.eccAdd(R,P,a,b,p)
        
        return R
    
    def eccMultiply_navie(self,k,P,a,b,p):
        ''' 穷算方法的椭圆曲线乘法： kP: k 次 P 相加
            a,b 为椭圆曲线参数， p 为模数
            
            复杂度为 O(k)
        '''
        R = P
        while k > 1:
            R = self.eccAdd(P,R,a,b,p)
            k -= 1
        return R
    def eccAttack(self,P,Q,a,b,p):
        ''' 求私钥 k,穷举攻击
        
        攻击者获得 P,Q，希望求得 k，Q = kP
        a,b 为椭圆曲线参数
        p 为模数
        
        计算时间复杂度 O(k)
    
        '''
        R = P
        k = 1
        while 1:
            if R == Q:
                break
            R = self.eccAdd(R,P,a,b,p)
            k += 1
            
        return k
    def isPrime(self,n):
        ''' 判断 n 是否为素数
        '''
        if n == 2: return True
        primeList = [2]
        for i in range(3,n+1):
            flag = True
            for j in range(0,len(primeList)):
                if i % primeList[j] == 0: 
                    flag = False
                    break
            if flag: primeList.append(i)
        #print(primeList)
        if primeList[-1] == n:
            return True
        return False
    def getPrimeList(self,a,b):
        ''' 获取 [a,b] 区间的所有素数
        '''
        if a > b: return []
        primeList = [2]
        pos = 0
        first = True
        if a == 2: first = False
        for i in range(3,b+1):
            flag = True
            for j in range(0,len(primeList)):
                if i % primeList[j] == 0:
                    flag = False
                    break
            if flag: 
                if i >= a and first: 
                    pos = len(primeList)
                    first = False
                primeList.append(i)
        
        return primeList[pos:]
            
        
    def millionarie(self,N):
        ''' 模拟姚氏百万富翁第一种解法
            N 是 素数的位数
        '''
        i , j = 3, 8
        print('millionarie A has ',i,'millioarie B has',j)
        x = int(math.pow(2,N))+random.randint(0,int(math.pow(2,N-1)))
        print('Bob choose x =',x)
        # alice public key e = 2, private key d = 4, module number n = 15
        e,d,n = 2,4,15
        k = self.quickReductionExp(x,e,n)
        print('k =',k,'Bob send k-j+1 =',k-j+1)
        y = [0]*10
        for u in range(0,10):
            y[u] = self.quickReductionExp(k-j+u,d,n)
        print('y1,y2,...y10 = ',y)
        
        primeList = self.getPrimeList(int(math.pow(2,N//2)),int(math.pow(2,N//2+1)))
        print('N/2 bit primeList: ', primeList)
        p = primeList[random.randint(0,len(primeList)-1)]
        print('Alice p=',p)  
        
        return
    
    def isQuadraticResidue(self,x,p):
        ''' 判定 x 是否为模 p 的二次剩余
        p 为素数
        x 为 p 以内的正整数
        '''
        flag = False
        for i in range(1,(p+1)//2):
            if i*i%p == x: 
                print('正根=',i,'负根=',p-i,'二次剩余为',x)
                flag = True
                
        return flag
    
    def showQR(self,p):
        '''' 显示 p 以内的所有平方剩余
        p 为素数
        容易证明，p 以内的平方剩余一共有 p//2 个
        '''
        for x in range(1,p):
            flag = False
            for i in range(1,(p+1)//2):
                if i*i%p == x: 
                    print(x,'为二次剩余','，正根=',i,'，负根=',p-i)
                    flag = True
            if flag == False: 
                print(x, '为平方非剩余')
                
        return
    def calcQR(self,p):
        ''' 计算 Zp* 内的平方剩余
        p 为素数
        返回值 r: 所有平方剩余形成的列表
        
        由平方剩余相关定理，保证了必然有 (p-1)/2 个平方剩余，而且就是前半部分的平方值取余数
        '''
        r = [0]*(p//2)
        for i in range(1,(p+1)//2):
            r[i-1] = i*i%p
            
        return r
    
    def linearCongruenceEquation(self,a,b,m):
        '''一次同余方程 ax=b(mod m), 求未知数 x
        返回值 r： 是一个列表，因为可能有多个解
        '''
        g = self.gcd(a,m)
        aa = int(a/g) # a'
        mm = int(m/g) # m'
        print(aa,mm)
        
        aaInv = self.getMulInverse(aa,mm) # a' 的逆元
        
        r = [0]*g
        r[0] = int(b/g)*aaInv%mm
        for i in range(1,g):
            r[i] = r[i-1]+mm
        
        return r
    
    def f(self,p):
        ''' 输出每个元素形成的循环子群
        元素1： 元素1 的循环子群
        '''
        for i in range(1,p):
            r = i
            print(r,':',end='{ ')
            while 1:
                print(r,end=' ')
                r = r*i%p
                if r==i: break
            
            print('}\n')
        
        return
    
    def isRoot(self, g,n):
        ''' 判定 g 是否为 模 n 的原根
        '''
        t = [g%n]
        i = 0
        while 1:
            t += [t[i]*g%n]
            i += 1
            print(t)
        
        
        
d = discreteMathematic_0110()
        
