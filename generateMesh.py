'''
points = []
results = []
validationPoints = []
validationResults = []
inputData = open('test.dta', 'r')
'''
f = open('Mesh2.node','w')
firstString = '400 2 2 0 \n'
f.write(firstString)
for i in range(400):
	vertexNumber = i + 1
	mod = i % 20
	width = -1.0 + 0.05 + 2 * (mod/20.0)
	quot = i / 20
	height = 1.0 - 0.05 - 2 * (quot/20.0)
	velx = -width
	vely = -height
	string = str(vertexNumber) + ' ' + str(width) + ' ' + str(height) + ' ' + str(velx) + ' ' + str(vely) + '\n' 
	f.write(string)
#counter = 0
'''
for line in inputData:
	#print line
	#counter += 1
	words = line.split()
	x1 = float(words[0])
	x2 = float(words[1])
	x3 = float(words[2])
	if x1 == 0:
		string = '+1 1:' + words[1] + ' 2:' + words[2] + '\n'
	else:
		string = '-1 1:' + words[1] + ' 2:' + words[2] + '\n'
	print words
	f.write(string) # python will convert \n to os.linesep

	
inputData.close()
'''
f.close() # you can omit in most cases as the destructor will call if