#/usr/bin/env awk
{
    count[NR] = $1;
}
END {
    if (NR % 2) {
        print count[(NR + 1) / 2];
    } else {
        print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0;
    }
}
