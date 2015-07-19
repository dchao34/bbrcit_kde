CREATE TEMPORARY VIEW E AS
SELECT 
    eid,
    rf_useopt_score,
    rf_dvsdstar_sigmc_score,
    bad,
    data_label
FROM (ONLY Event INNER JOIN EventMetaData USING (eid));

CREATE TEMPORARY VIEW C AS
SELECT *
FROM ONLY Candidate INNER JOIN OptCandidateIdx USING (eid, idx);

CREATE TEMPORARY VIEW F AS
SELECT 
    rf_useopt_score,
    rf_dvsdstar_sigmc_score,
    tag_lp3
FROM E INNER JOIN C USING (eid)
WHERE
    bad=0 AND
    data_label=0
;

\copy (SELECT * FROM F) TO 'fitdata_offpeak.csv' WITH CSV HEADER;
