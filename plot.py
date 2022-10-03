import matplotlib.pyplot as plt
import os

if os.path.exists('./phys_output.dat'):
  
  time = []
  etot = []
  epot = []
  ekin = []


  file_output = open('phys_output.dat', 'r')
  next(file_output)
  next(file_output)
  for line in file_output:
    values = [float(s) for s in line.split()]
    time.append(values[0])
    etot.append(values[1])
    epot.append(values[2])
    ekin.append(values[3])

  plt.plot(time, etot, '-k', label = '<E$_{tot}$>',c ='black')
  plt.plot(time, epot, '--k', label = '<E$_{pot}$>', c ='blue')
  plt.plot(time, ekin, '-.k', label = '<E$_{kin}$>', c = 'red')
    
  plt.xlabel('$t_{step}$', fontsize = 12)
  plt.ylabel('V$_{ij}$', fontsize = 12)
    
  plt.title('Mean Energy', fontsize = 20)
  plt.legend()
  plt.savefig('mean_energy.png')
  plt.show()

