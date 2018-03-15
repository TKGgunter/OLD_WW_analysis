datatime=`date +%Y_%m_%d`
tar --exclude='Analysis_13TeV/data' --exclude='Analysis_13TeV/plots' -cvf ~/STORAGE/WW_analysis_$datatime.tar Analysis_13TeV
