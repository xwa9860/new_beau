lines = open('fuel_mcnp', 'r+').readlines()
fuel_eq = []
for line in lines:
	fuel_eq.append(line.split()[0].split('.')[0] + ' %E\n' %(float(line.split()[1]) + float(line.split()[3])))
output = open('fuel_eq', 'w+')
output.write(''.join(fuel_eq))
