import re
from lib import prep

uniprots = []
with open('uni.fasta', "r") as f:
        content = f.read()
        identifiers = re.findall(r"sp\|([^|]+)\|", content)
        for identifier in identifiers:
            uniprots.append(identifier)

with open('important.txt', "w") as f:
    for i in uniprots:
        a = prep.query_uniprot_for_glycosylation_locations(i)['glycosylations']
        if a!= []:
            f.write(i+"\n")
            print(i)
