.headers on
.mode csv
select t-s as stem, s, t, count(*) as dim from SteenrodMRes_generators GROUP BY s, t ORDER BY stem;