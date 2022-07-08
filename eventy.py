L2=140
L3=300
LHL=3000
brHWWbb=.21*.58
xsec = .266
ib=input("Is this L2, L3, or LHL? ")
while ib != "no":

	ia=input("Branching ratio: ")
	ia=float(ia)


	if ib == "L2":
		eventy = xsec*L2*ia*brHWWbb
		print(eventy)
	if ib == "L3":
		eventy = xsec*L3*ia*brHWWbb
		print(eventy)
	if ib == "LHL":
		eventy = xsec*LHL*ia*brHWWbb
		print(eventy)
	ib=input("Is this L2, L3, or LHL, type no to stop? ")
