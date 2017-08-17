import xlrd

def read_xlrd(filename):
	workbook = xlrd.open_workbook(filename)
	booksheet = workbook.sheet_by_name('Sheet1')
	nameMap = {}
	for row in range(booksheet.nrows):
		xlsx_item = []; dic={}
		for col in range(booksheet.ncols):
			cel = booksheet.cell(row,col)
			val = cel.value
			xlsx_item.append(val.strip().split())	
		for item in xlsx_item[1]:
			dic[item] = xlsx_item[0][0]
		nameMap.update(dic)
	for item in nameMap.items():
		print(item)
	return nameMap
	
if __name__ == '__main__':
	read_xlrd('test.xlsx')