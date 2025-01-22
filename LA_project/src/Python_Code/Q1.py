import numpy as np
import pandas as pd

class Matrix_Operation:
    def __init__(self,matrix):
        self.__matrix=np.array(matrix)

    def set_matrix(self,array):
        array=np.array(array)
        self.__matrix=array

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






    








a=np.array([[1,0,0],[0,1,0],[0,0,1]])
Matrix_op=Matrix_Operation(a)
# Matrix_op.show_matrix()
# print(Matrix_op.transpose())
Matrix_op.inverse()
Matrix_op.show_matrix()







        
            
    




