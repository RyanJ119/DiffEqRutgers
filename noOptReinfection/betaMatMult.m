mat = [1 2 3; 4 5 6; 7 8 9]
mat2 = [1 2 1]'
vector = [1; 5; 7] %susceptibles



result2 = mat.*vector.*vector

result3 = sum((mat.*vector.*vector), 2)



result*mat2'