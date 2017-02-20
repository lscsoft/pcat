#!/bin/bash  

for i in "$@"
do
case $i in
    -d=*|--date=*)
    DATE="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

source ~detchar/opt/gwpysoft/bin/activate

if [ ! -n "$DATE" ]
then
  #today=`date --utc +%Y-%m-%d`
  day=`date --utc --date='yesterday 00:00' +%Y-%m-%d`
else
  day=$DATE
fi

day_gps=`lalapps_tconvert $day`
final_time=`expr $day_gps + 86400`
initial_time=$day_gps

pcat-multi --IFO L --start=$initial_time --end=$final_time -t PCAT_time.config -f PCAT_frequency.config -n $day

exit
