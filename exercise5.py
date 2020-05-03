seq = getSequence('MN908947', 'genbank', DNA_Alphabet)

#So this is what is causing alot of trouble right now....
print(seq.sequence)


n_a=seq.count("A")
n_c=seq.count("C")
n_g=seq.count("G")
n_t=seq.count("T")
total=n_a+n_c+n_g+n_t

per_a=round(n_a*100/total)
per_c=round(n_c*100/total)
per_g=round(n_g*100/total)
per_t=round(n_t*100/total)
print(per_a)
print(per_c)
print(per_g)
print(per_t)




seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(0,True)
x1=trans
print(trans)


seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(0,False)
x2=trans
print(trans)


seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(1,True)
x3=trans
print(trans)


seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(0,False)
x4=trans
print(trans)


seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(1,False)
x5=trans
print(trans)


seq=Sequence(seq,DNA_Alphabet)
# print(seq)
seq.translateDNA(2,False)
x6=trans
print(trans)



start_1=[i for i,v in enumerate(x1) if v=="M"]
# print(start_1)
end_1=[i for i,v in enumerate(x1) if v=="*"]
# print(end_1)



# orf=0
# 开始位置列表取第i个元素时，
# 结束位置列表找到第一个大于开始位置列表取第i个元素
# 如果结束位置列表第一个大于开始位置列表取第i个元素的元素-开始位置列表第i个元素>100
# orf+1