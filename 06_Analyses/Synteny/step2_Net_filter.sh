if [ "$#" -lt 2 ]; then
    echo "Usage: sh step2_Net_filter.sh ../UCSC.target.filter.net ../../query.sizes "
    exit 1
fi

net=$1
size=$2

perl netFill.pl $net > net.lst
perl chgOri.pl $size net.lst > net.chg.lst
awk '($8>500000 || $9>500000)' net.chg.lst > net.lst.filter
perl -ne '@t=split /\s+/;pop @t;pop @t;print join("\t",@t),"\n"' net.lst.filter   | grep -vi Scaffold | grep -v NW_  > net.lst.filter2
echo -e "s1chr\tStart_1\tEnd_1\ts2chr\tStart_2\tEnd_2" > header
cat header net.lst.filter2 > input.lst
perl -ne '@t=split /\s+/;pop @t;pop @t;print join("\t",@t),"\n"' net.chg.lst | grep -vi Scaffold > net.lst.filter2.all
cat header net.lst.filter2.all > input.lst.all
