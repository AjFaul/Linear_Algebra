import numpy as np
import pandas as pd
import math as m


class Matrix_Operation:
    def __init__(self,matrix):
        self.__matrix=np.array(matrix)

    def set_matrix(self,array):
        array=np.array(array)
        self.__matrix=array
    
    def get_matrix(self):
        return self.__matrix

    def show_matrix(self,matrix=None):
        if matrix is None :
            print(matrix)
        print(self.__matrix)

    def transpose(self):
        (rows,cols)=(self.__matrix.shape)
        transpose_arr=np.zeros((cols,rows))
        for i in range(rows) :
            for j in range(cols):
                transpose_arr[j][i]=self.__matrix[i][j]
        self.__matrix=transpose_arr
        return transpose_arr
    
    def determinant(self, array=None):
        if array is None:
            array = self.__matrix
            
        n = len(array)
        if n != array.shape[1]:
            raise ValueError("Matrix must be square to compute determinant.")
            
        if n == 1:
            return array[0, 0]
            
        if n == 2:
            return array[0, 0] * array[1, 1] - array[0, 1] * array[1, 0]
            
        det = 0
        for col in range(n):
            sub_matrix = np.delete(np.delete(array, 0, axis=0), col, axis=1)
            det += ((-1) ** col) * array[0, col] * self.determinant(sub_matrix)
        return det
    
    def adjugate(self):
        n = len(self.__matrix)
        if n != self.__matrix.shape[1]:
            raise ValueError("Matrix must be square to compute adjugate.")
        adj = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                sub_matrix = np.delete(np.delete(self.__matrix, i, axis=0), j, axis=1)
                adj[i][j] = ((-1)**(i + j)) * self.determinant(sub_matrix)
        self.__matrix=adj
        return self.transpose()
    
    def inverse(self):
        det = self.determinant()
        if det == 0:
            raise ValueError("Matrix is singular and cannot be inverted.")
        adj = self.adjugate()
        inv = adj / det
        self.__matrix=inv
        return inv
    
    def multiply_matrices(self, array1, array2):
        array1 = np.array(array1)
        array2 = np.array(array2)
            
        if array1.shape[1] != array2.shape[0]:
            raise ValueError("Incompatible dimensions for matrix multiplication.")
            
        rows_array1 = array1.shape[0]
        cols_array2 = array2.shape[1]
        cols_array1 = array1.shape[1]
            
        result = np.zeros((rows_array1, cols_array2))
            
        for i in range(rows_array1):
            for j in range(cols_array2):
                sum = 0
                for k in range(cols_array1):
                    sum += array1[i][k] * array2[k][j]
                result[i][j] = sum
            
        return result

    def A_back_slash_b(self, b):
        a=self.__matrix
        a_t=self.transpose()
        self.__matrix=self.multiply_matrices(a_t,a)
        at_a_inv=self.inverse()
        a_cross=self.multiply_matrices(at_a_inv,a_t)
        return self.multiply_matrices(a_cross,b)
    

class Jacobi:
    
    def rotation_matrices(self,size, theta,p,q):
        # print(f"theta is {theta}")
        E=np.eye(size)
        E[p][q]=m.sin(theta)
        E[p][p]=m.cos(theta)
        E[q][q]=m.cos(theta)
        E[q][p]=-1 * m.sin(theta)
        # print(E)
        return E 

    def max_abs_element(self,matrix,num_row,num_col):
        max_num=[-1000000,0,0]
        for i in range(num_row):
            for j in range(num_col):
                if max_num[0]<abs(matrix[i][j]) and i!=j:
                    max_num[0]=abs(matrix[i][j])
                    max_num[1]=i
                    max_num[2]=j
        return max_num
  
    def find_theta(self,matrix):
        j=Jacobi()
        shape=list(matrix.shape)
        max_element=j.max_abs_element(matrix,shape[0],shape[1])
        p=max_element[1]
        q=max_element[2]
        if (matrix[p][p] - matrix[q][q]) ==  0 :
            theta= m.pi * 1/4 * (matrix[p][q] / abs(matrix[q][p]))
        else:
            theta=(1/2) * m.atan(-1 * 2 * matrix[p][q] / (matrix[p][p] - matrix[q][q] ))
        return theta

    def create_E(self,matrix):
        shape=list(matrix.shape)
        theta=self.find_theta(matrix)
        max_num=self.max_abs_element(matrix,shape[0],shape[1])
        # print(f"Hi p is {max_num[1]} and q is {max_num[2]}")
        return self.rotation_matrices(shape[1],theta,max_num[1],max_num[2])

        
def jacobi_method(matrix,num_of_iteration):
    j=Jacobi()
    list_of_E=[]
    m=Matrix_Operation(matrix)
    for i in range(num_of_iteration):
        e=j.create_E(m.get_matrix())
        E=Matrix_Operation(e)
        list_of_E.append(E)
        m.set_matrix(m.multiply_matrices(E.transpose(),m.get_matrix()))
        m.set_matrix(m.multiply_matrices(m.get_matrix(),E.transpose()))
    return [list_of_E,m]
        
        

    
    
    
    
    
    
    
    
    
    
    

A=np.array([[1,7,3,5],[7,4,6,2],[3,6,0,2],[5,2,2,-1]])
ans=jacobi_method(A,100)
ans[1].show_matrix(ans[1].get_matrix())



    
    
    
    