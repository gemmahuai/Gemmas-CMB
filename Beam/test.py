import numpy as np
import csv

amp = np.array([0.1, 1.0, 3.0, 5.0, 10.0, 20.0])
k_in = [1, 10, 20, 50]
error = np.zeros((len(k_in),len(amp)))
error = np.array([[1,2,3,4,5,6],[2,3,4,5,6,7],[6,6,6,6,6,6], [10,11, 12, 13, 14, 15]])
print('k_in = ', np.repeat(k_in, len(amp)))
print('amp = ', np.tile(amp, len(k_in)))
print('error = ', (error**2).flatten())
#print(screen)

### write to a csv file
# a = amp #amplitude of errors
# b = screen #leakag
# c = np.array([10,50])
# print(zip(a,b,c))
# with open('data_test.csv', 'w') as f:
#    writer = csv.writer(f, delimiter='\t')
#    writer.writerows(zip(a,b,c))