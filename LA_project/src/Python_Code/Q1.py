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

    def A_back_slash_b(self, b):
        a=self.__matrix
        a_t=self.transpose()
        self.__matrix=self.multiply_matrices(a_t,a)
        at_a_inv=self.inverse()
        a_cross=self.multiply_matrices(at_a_inv,a_t)
        return self.multiply_matrices(a_cross,b)

    def LU_decomposition(self):
        n,m =self.__matrix.shape
        if n is not m:
            raise ValueError("Is not Square")
        L=np.eye(n)
        U=self.__matrix.astype(float)
        for k in range(n-1):
            if U[k,k]==0:
                raise ValueError("is singular")
            for i in range(k+1,n):
                L[i,k]=U[i,k]/U[k,k]
                U[i,:]-=L[i,k]*U[k,:]
        return L,U
    
    def solve_lu(self,b):
        L,U=self.LU_decomposition()
        n=len(b)
        y=np.zeros(n)
        for i in range(n):
            y[i]=b[i]- np.dot(L[i,:i],y[:i])
        
        n=len(y)
        x=np.zeros(n)
        for i in range(n-1,-1,-1):
            x[i]= (y[i]- np.dot(U[i,i+1:] , x[i+1:]))/U[i,i]
        return x
        



A=np.array([[1,7,3,5],[7,4,6,2],[3,6,0,2],[5,2,2,-1]])
b=np.array([[1],[2],[5],[-1]])
#init matrix
Matrix_op=Matrix_Operation(A)

# A/b
print(Matrix_op.A_back_slash_b(b))

# x=Matrix_op.solve_lu(b)
# print(x)









        
            
    




