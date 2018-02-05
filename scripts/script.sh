declare -a list=($(find ~/Dropbox/PapersOR/PM/implementation/resources/instances/ -regex '.+[01]+[045]0\/.+[16]\.dat'))
for m in   2 4; do
for alpha in (0.7 0.75 0.8 0.85 0.9)  ; do
for file in $list;do
./PM  -S 1 -a $alpha -f 5  $file $m;
echo $m;
done
done
done
