--psql bfw

set search_path TO public, owi;

\copy (select peri, rw, hw, pbfl, tlfl, azi, dist, bart, bhd, d03h, hoehe, kronho, nutzart, ea, ba from bam where peri in (6,7)) to /tmp/bam.csv csv header
\copy (select peri, rw, hw, pbfl, datum, kg_wald, seehoehe, rueck, vorrueck from prf where peri > 2) to /tmp/prf.csv csv header
\copy (select peri, rw, hw, pbfl, tlfl, ba, ea, hangneig, kg_wald, neigricht, relief, wugeb, bogrup, vegtyp  from tfl where peri > 2) to /tmp/tfl.csv csv header


--mv /tmp/bam.csv ./dat
--mv /tmp/prf.csv ./dat
--mv /tmp/tfl.csv ./dat
--bzip2 ./dat/bam.csv ./dat/prf.csv ./dat/tfl.csv
