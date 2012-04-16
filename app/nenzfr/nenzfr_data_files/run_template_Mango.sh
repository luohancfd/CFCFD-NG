# This file should not need to be modified.
# Luke Doherty
# 16-04-2012


nohup nenzfr.py --gas=$gasName --T1=$T1 --p1=$p1 --Vs=$Vs --pe=$pe --chem=$chemModel --area=$areaRatio --job=$jobName --cfile=$contourFileName --gfile=$gridFileName --exitfile=$exitSliceFileName $blockMarching --nni=$nni --nnj=$nnj --nbi=$nbi --nbj=$nbj --bx=$bx --by=$by --max-time=$max_time --max-step=$max_step --Twall=$Tw --BLTrans=$BLTrans --TurbVisRatio=$TurbVisRatio --TurbIntensity=$TurbInten --CoreRadiusFraction=$CoreRadiusFraction > LOGFILE &


