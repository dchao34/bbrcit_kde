Data source
-----------

1. Download data from database using `download_generic.sql`. 
   + Output: `fitdata_generic.csv`. 

2. Partition training/testing set. 
   + Program: `split_generic.cc`
   + Input: `fitdata_generic.csv`. 
   + Output: `train_generic.csv` and `test_generic.csv`. 

3. Perform event weighting by down sampling. 
   + Program: `weighted_down_sample.cc`. 
   + Input: `train_generic.csv` and `test_generic.csv`. 
   + Output: `train.csv` and `test.csv`. 

4. Extract columns and event types that are relevant to KDE. 
   + Program: `prepare_kde_data.cc`. 
   + Input: `train.csv` or `test.csv`, and event types. 
   + Output: Only the following types of inputs are defained.
     1. If event types are 1, 2, 3, 4, *and* 5 then the results are 
        `train.kde.csv` or `test.kde.csv`. 
     2. If event types are 1, 2, 3, 4, *or* 5, *and* if `train.csv` 
        is given, then the result is `evttype{t}.csv`, where `t` is
        the specified event type. 

5. Subsample data for KDE cross validation. 
   + Program: `subsample.cc`. 
   + Input: `evttype{t}.csv`, `t`=1, 2, 3, 4, or 5. 
   + Output: `evttype{t}.cv.csv`. 
