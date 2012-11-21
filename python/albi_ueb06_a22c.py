
query = open('query.fa', 'r')
db = open('library.fa','r')

# lese query sequenz
query_string = query.readlines()[1][:-1]

count = 0

for line in db:
    if( query_string in line):
        ++count
    

print (count)

query.close()
db.close()
        
