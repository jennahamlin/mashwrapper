import os


#dirpath = os.path.dirname(os.path.realpath(__file__))
#filepath = os.path.join(dirpath, 'test-data/file1.txt')
#print(dirpath)
#print(filepath)
#open(filepath)  # or whatever you use to open the file
#with open(filepath) as f:
#    contents = f.read()
#    print(contents)

dirpath = os.path.dirname(os.path.realpath(__file__))
filepath1 = os.path.join(dirpath, 'test-data/subERR125190_1.fastq.gz')
filepath2 = os.path.join(dirpath,'test-data/subERR125190_2.fastq.gz')
filepath3 = os.path.join(dirpath, 'bin/test-data/myMashDatabase.msh')

## open(filepath)  # or whatever you use to open the file
with open(filepath3) as myDB:
    contents = myDB.read()
    print(contents)
