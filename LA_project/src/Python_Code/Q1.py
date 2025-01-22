import numpy as np
import pandas as pd

class Matrix_Operation:
    def __init__(self,matrix):
        self.__matrix=np.array(matrix)

    def show_matrix(self):
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
    


    








a=np.array([[1,2,3],[0,2,2],[0,0,7]])
Matrix_op=Matrix_Operation(a)
Matrix_op.show_matrix()
print("\n------------------\n")
print(f"deteminan is {Matrix_op.determinant()}")






        
            
    




