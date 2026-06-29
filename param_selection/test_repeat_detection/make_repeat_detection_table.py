#!/usr/bin/env python3

import os
import sys

import polars as pl


input_fpath = os.path.abspath(sys.argv[1])
out_fpath = os.path.abspath(sys.argv[2])
finder = sys.argv[3]


# >>> Proceed >>>

df = pl.read_csv(input_fpath, separator='\t')
print('INPUT TABLE:')
print(df.head())

df = df.with_columns(
    (
        'MR_l' + pl.col('repeat_len').cast(pl.String) \
        + '_' \
        + 'r' + pl.col('replicate_idx').cast(pl.String) \
        + '_' \
        + 't' + pl.col('rate').cast(pl.String) \
    ).alias('repeat_id')
)

if finder == 'phraider':
    print('FINDER is `phraider`')
    print('Removing reverse-complement repeats...')
    df = df.filter(
        ~(
            (pl.col('value_type') == 'true') & (pl.col('strand') == '-')
        )
    )
# end if

true_df = df.filter(
    pl.col('value_type') == 'true'
)
pred_df = df.filter(
    pl.col('value_type') == 'pred'
)

all_repeat_ids = frozenset(
    true_df['repeat_id']
)

chr_list = list()
repeat_id_list = list()
true_start_coord_list = list()
true_end_coord_list = list()
pred_start_coord_list = list()
pred_end_coord_list = list()
repeat_found_list = list()

for repeat_id in all_repeat_ids:
    curr_true_df = true_df.filter(pl.col('repeat_id') == repeat_id)
    curr_pred_df = pred_df.filter(pl.col('repeat_id') == repeat_id)

    for true_row in curr_true_df.to_dicts():

        curr_pred_detect_df = curr_pred_df.filter(
            (pl.col('repeat_id') == repeat_id) \
          & (pl.col('chr') == true_row['chr']) \
          & (pl.col('end_coord') >= true_row['start_coord']) \
          & (pl.col('start_coord') <= true_row['end_coord'])
        )

        if curr_pred_detect_df.height > 0:
            found_value = 1
            for pred_row in curr_pred_detect_df.to_dicts():
                chr_list.append(true_row['chr'])
                repeat_id_list.append(repeat_id)
                true_start_coord_list.append(true_row['start_coord'])
                true_end_coord_list.append(true_row['end_coord'])
                pred_start_coord_list.append(pred_row['start_coord'])
                pred_end_coord_list.append(pred_row['end_coord'])
                repeat_found_list.append(found_value)
            # end for
        else:
            found_value = 0
            chr_list.append(true_row['chr'])
            repeat_id_list.append(repeat_id)
            true_start_coord_list.append(true_row['start_coord'])
            true_end_coord_list.append(true_row['end_coord'])
            pred_start_coord_list.append(None)
            pred_end_coord_list.append(None)
            repeat_found_list.append(found_value)
        # end if
    # end for
# end for

repeat_detect_df = pl.DataFrame({
    'chr': chr_list,
    'repeat_id': repeat_id_list,
    'true_start_coord': true_start_coord_list,
    'true_end_coord': true_end_coord_list,
    'pred_start_coord': pred_start_coord_list,
    'pred_end_coord': pred_end_coord_list,
    'repeat_detected': repeat_found_list,
})

repeat_detect_df = repeat_detect_df.join(
    true_df,
    on='repeat_id',
    how='inner'
).select(
    pl.col(
        'chr',
        'repeat_id',
        'true_start_coord',
        'true_end_coord',
        'pred_start_coord',
        'pred_end_coord',
        'repeat_detected',
        'repeat_len',
        'replicate_idx',
        'rate'
    )
).unique()

print('OUTPUT TABLE:')
print(repeat_detect_df.shape)
print(repeat_detect_df.head())

repeat_detect_df.write_csv(
    out_fpath,
    separator='\t',
    null_value='NA'
)

print(out_fpath)
sys.exit(0)
