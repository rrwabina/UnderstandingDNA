import sys
import string
import re
#import itertools
import nltk

class ConvertCNF():
	def __init__(self, rules, startSymbol):
		self.rules = rules
		self.startSymbol = startSymbol
		self.terminal = {}
		self.temp_list = {}
		self.eliminateEpsilon()
		self.eliminateVariableUnit()
		self.moveTerminalToUnit()
		self.replaceLongProd()
		self.startSymbolAdd()

	def startSymbolAdd(self):
		for item in self.rules[self.startSymbol]:
			if self.startSymbol in item:
				prev_start = self.startSymbol
				while self.startSymbol in self.rules:
					self.startSymbol = self.startSymbol+'0'
				self.rules[self.startSymbol] = self.rules[prev_start]
			else:
				pass

	def eliminateEpsilon(self):
		for key,value in self.rules.items():
			if 'ε' in value:
				for key2, value2 in self.rules.items():
					if key in str(value2):
						if key in value2:
							#self.rules[key2].append('ε')
							self.rules[key2] = self.rules[key2] + self.rules[key]
							self.rules[key2].remove('ε')
						else:
							for item in value2:
								new_prod = self.createProdCombinations(item, key, item.count(key))
								self.rules[key2] = self.rules[key2] + new_prod
				self.rules[key] = list(set(self.rules[key]))
				self.rules[key].remove('ε')

		for key, value in self.rules.items():
			for i in range(len(value)):
				self.rules[key][i] = self.rules[key][i].strip()

		for key, value in self.rules.items():
			value = set(value)
			value = list(value)
			self.rules[key] = value
		#print("After Epsilon Elimination: ")
		#print(self.rules)

	def createProdCombinations(self, item, nonterminal, count):
	        numset = 1 << count
	        new_prods = []

	        for i in range(numset):
	            nth_nt = 0
	            new_prod = ''
	            for s in item:
	            	if s == nonterminal:
	            		if i & (1 << nth_nt):
	            			new_prod = new_prod+s
	            		nth_nt += 1
	            	else:
	            		new_prod = new_prod+s
	            new_prods.append(new_prod)
	        return new_prods

	def eliminateVariableUnit(self):
		for key, value in self.rules.copy().items():
			flag=0
			for key2, value2 in self.rules.copy().items():
				if key in str(value2):
					if key in value2 and key==key2:
						self.rules[key].remove(key)
					if key in value2:
						self.rules[key2] = self.rules[key2] + self.rules[key]
						self.rules[key2].remove(key)
					for c in value2:
						if key in c:
							if key==c and flag!=2:
								flag=1
							else:
								flag=2
								break
					#print(flag)
			if flag==1 and key!=self.startSymbol:
				self.rules.pop(key, None)


		#print("After eliminate variable unit")
		#print(self.rules)

	def replaceLongProd(self):
		for key, value in self.rules.items():
			for i in range(len(value)):
				self.rules[key][i] = re.sub(' +',' ',value[i])

		for key, value in self.rules.copy().items():
			for key2, value2 in self.rules.copy().items():
				if key!=key2:
					for i in range(len(value)):
						for j in range(len(value2)):
							if value2[j] in value[i] and len(value2)==1:
								#print(value2)
								if len(value2[j].split(" "))!=2:
									self.rules[key][i]=self.rules[key][i].replace(value2[j], key2)
								#pass

		for key, value in self.rules.copy().items():
			for i in range(len(value)):
				if value[i] in self.rules:
					if all( tkey.startswith("'") for tkey in self.rules[value[i]]):
						self.rules[key][i] = self.rules[key][i].replace(value[i],self.rules[value[i]][0] )


		for key, value in self.rules.copy().items():
			for i in range(len(value)):
				splitted = value[i].split(" ")
				while len(self.rules[key][i].split(" "))>2:
					item = self.getNewNTSymbol()
					if self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1] not in self.temp_list:
						self.temp_list[self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1]] = item
					
					self.rules[self.temp_list[self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1]]] = [self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1]]
					self.rules[key][i] = self.rules[key][i].replace(self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1], self.temp_list[self.rules[key][i].split(" ")[0]+" "+self.rules[key][i].split(" ")[1]])

	def getNewNTSymbol(self):
		nt = list(string.ascii_uppercase)
		new_nt = ''
		for item in nt:
			if item in self.rules:
				continue
			else:
				new_nt = item
				break
		return new_nt

	def moveTerminalToUnit(self):
		for key, value in self.rules.items():
			for i in range(len(value)):
				if re.compile(r'\'+\w+\'').search(value[i]) and any(c in value[i] for c in self.rules):
					match = re.findall(r'\'+\w+\'', value[i])
					for i in match:
						self.terminal[i] = key
		for item, value in self.terminal.items():
			dict_new_nt = item.upper()
			dict_new_nt = dict_new_nt.replace("'","")
			while dict_new_nt in self.rules:
				dict_new_nt = self.getNewNTSymbol()
			self.rules[dict_new_nt] = [item]
			self.terminal[item] = dict_new_nt
		#print("After move terminal to unit:")
		#print(self.rules)

class CYK():
	def __init__(self, path, string_path):
		self.temp_table = []
		self.path = path
		self.string_path = string_path
		self.startSymbol = ''
		self.rules = {}
		self.words = []
		print("\nInput Grammar:")
		self.readGrammar()
		self.readString()
		self.converted= ConvertCNF(self.rules, self.startSymbol)
		i=1
		for item in self.words:
			self.word = item
			print("\n"+str(i)+"----------------------------")
			print("\nFinal Grammar:")
			i+=1
			self.readOutput()
			print(self.word)
			self.parser()

	def readGrammar(self):
		f = open(self.path)
		self.startSymbol = f.readline().rstrip()
		for content in f:
			print(content.strip("\n"))
			content = content.rstrip()
			rule = content.split(" -> ")

			if rule[0] not in self.rules:
				self.rules[rule[0]] = rule[1].split(" | ")
			else:
				self.rules[rule[0]] += (rule[1].split(" | "))

	def readString(self):
		f = open(self.string_path)
		for content in f:
			content = content.rstrip()
			self.words.append(content)

	def readOutput(self):
		str1=''
		for key, value in self.converted.rules.items():
			str1 += str(key) + " -> "+ ' | '.join(map(str, value))+'\n'
			#print( key + " -> "+ ' | '.join(map(str, value))
		print(str1)

	
		groucho_grammar = nltk.CFG.fromstring(str1)
		sent = self.word.split(" ")
		parser = nltk.ChartParser(groucho_grammar)
		print("\nParse Tree:")
		for tree in parser.parse(sent):
			print(tree)
		

	def parser(self):
		self.copy_table = []
		wordList = self.word.split(" ")
		length = len(wordList)
		#print(length)
		self.parse_table = [[[] for x in range(length - y)] for y in range(length)]

		for i, word in enumerate(wordList):
			for key, value in self.converted.rules.items():
				for item in value:
					if item== "'"+word+"'":
						self.parse_table[0][i].append(key)
						self.copy_table.append((key, item))

		i=0
		for possible_words in range(2, length + 1):
			for current_cell in range(0, length - possible_words + 1):
				for left_side in range(1, possible_words):
					right_side = possible_words - left_side

					left_cell = self.parse_table[left_side - 1][current_cell]
					right_cell = self.parse_table[right_side - 1][current_cell + left_side]

					left_flag=0
					right_flag=0
					for key, value in self.converted.rules.items():
						for item in value:
							subitem = item.split(" ")
							for left in left_cell:
								if left==subitem[0]:
									left_flag=1
									for right in right_cell:
										if right==subitem[1]:
											right_flag=1
											if left_flag==1 and right_flag==1:
												print(key)
												self.parse_table[possible_words - 1][current_cell].append(key)
												self.copy_table.append((key, [subitem[0], subitem[1]], (left_side - 1, current_cell), (right_side - 1, current_cell + left_side), (possible_words - 1, current_cell) ))
										left_flag=0
										right_flag=0
		#self.copy_table = self.parse_table

		for item in self.copy_table:
			if item[0] != self.startSymbol:
				if all(item[0] not in item2[1] for item2 in self.copy_table):
					self.copy_table.remove(item)
		
		print("\nCKY Chart")
		for item in reversed(self.parse_table):
			print(' __ '.join(map(str, item)))
		print(" ____ ".join(map(str, self.word.split(" "))))
		
		
		if self.parse_table[-1][0]:
			if self.startSymbol in self.parse_table[-1][0]:
				#print("\n")
				self.copy_table = self.copy_table[::-1]
				#print("\n")
				self.copy_table_reversed = self.copy_table[::-1]
				print("----")
				print(self.copy_table_reversed)
				print("----")
				a = self.generate_tree(self.copy_table[0])
				print("Parse Tree:")
				print(a)
				#for item in self.copy_table:
				#	print(item)
				print("\nAccepted String")
		else:
			print("\nRejected String")

	def generate_tree(self, node):
		#print(node)
		if isinstance(self.find(self.copy_table, node)[1], (list,)):
			item = self.find(self.copy_table, node)[1]
			item_x = self.find(self.copy_table, node)[2]
			item_y = self.find(self.copy_table, node)[3]

			if item_x[0]==0:
				x=self.copy_table_reversed[item_x[1]]
			else:
				x=self.find_tuple(self.copy_table, item_x)

			if item_y[0]==0:
				y=self.copy_table_reversed[item_y[1]]
			else:
				y=self.find_tuple(self.copy_table, item_y)

			return (node[0], self.generate_tree(x), self.generate_tree(y))
		else:
			return (node[0],  self.find(self.copy_table, node)[1])

	def find_tuple(self, list_of_tuples, value):
		val = []
		for x in list_of_tuples:
			if x[4]==value:
				val = x
				break
		return val

	def find(self, list_of_tuples, value):
		val = []
		for x in list_of_tuples:
			if x==value:
				val = x
				break
		return val

if __name__ == "__main__":
    CYK(sys.argv[1], sys.argv[2])
