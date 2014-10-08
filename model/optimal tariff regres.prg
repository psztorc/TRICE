matrix(10,10) results
smpl @all
!m=1
for !n = 1 to 7
ls t{!n}  w{!n}  c
'ls t{!n}  w{!n} w{!n}^2   w{!n}^3 c
results(!m,1)=@coefs(1)
'results(!m,2)=@coefs(2)
'results(!m,3)=@coefs(3)
'results(!m,2)=@tstat(2)
'results(!m,3)=@r2




!m=!m+1
next
show results
