import matplotlib.pyplot as plt
import os

if os.path.exists('./phys_output.dat'):
  
  time = []
  etot = []
  epot = []
  ekin = []
    
  for line in open('phys_output.dat', 'r'):
    values = [float(s) for s in line.split()]
    time.append(values[0])
    etot.append(values[1])
    epot.append(values[2])
    ekin.append(values[3])

  plt.plot(time, etot, '-k', label = '<etot>')
  plt.plot(time, epot, '--k', label = '<epot>')
  plt.plot(time, ekin, '-.k', label = '<ekin>')
    
  plt.xlabel('time', fontsize = 12)
  plt.ylabel('', fontsize = 12)
    
  plt.title('Mean Energy', fontsize = 20)
  plt.legend()
  plt.savefig('mean_energy.png')
  plt.show()

