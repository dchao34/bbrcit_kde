CREATE TEMPORARY VIEW E AS
SELECT 
    eid,
    rf_useopt_score,
    rf_dvsdstar_sigmc_score,
    event_weight,
    mc_evttype,
    bad,
    data_source,
    ml_sample
FROM (McEvent INNER JOIN EventMetaData USING (eid));

CREATE TEMPORARY VIEW C AS
SELECT *
FROM McCandidate INNER JOIN OptCandidateIdx USING (eid, idx);

CREATE TEMPORARY VIEW F AS
SELECT 
    rf_useopt_score,
    rf_dvsdstar_sigmc_score,
    tag_lp3,
    event_weight,
    mc_evttype
FROM E INNER JOIN C USING (eid)
WHERE
    bad=0 AND
    data_source>1 AND 
    ml_sample=-1
;

\copy (SELECT * FROM F) TO 'fitdata.csv' WITH CSV HEADER;
