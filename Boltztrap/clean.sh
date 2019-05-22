dn="dn"
up="up"
for name in K*
do
cd $name
rm -r $name$dn $name$up
cd ..
done
