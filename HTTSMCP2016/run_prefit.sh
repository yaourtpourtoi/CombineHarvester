# for ch in tt mt et em
#     do python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_${ch}_2016_1_ -f output/301118_newbbb/shapes.root --ratio --x_axis_min 50 --combined_yrs
# done

# for era in 2016 2017
#     do for ch in tt mt et em
#         do for bin in {2..6}
#             do python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_${ch}_${era}_${bin}_ -f shapes.root -f_alt shapes_ps.root --ratio  --log_y --proper_errors_asym 
#         done
#     done
# done

for era in 2016 2017
    do for ch in tt mt et em
        do python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_${ch}_${era}_1_ -f shapes.root --ratio --x_axis_min 50 
    done
done

for era in 2016 2017
    do for ch in tt mt et em
        do for bin in {2..6}
            do python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_${ch}_${era}_${bin}_ -f shapes.root -f_alt shapes_ps.root --ratio --log_y
        done
    done
done

