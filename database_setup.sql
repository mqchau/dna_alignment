drop schema public cascade;
create schema public;

CREATE TABLE IF NOT EXISTS reference_hash (
    mer     VARCHAR(20),
    location INTEGER[],
    PRIMARY KEY (mer)
);

CREATE INDEX mer_hash ON reference_hash USING hash (mer);


-- This table is used to store the result after we do local alignment after hashing match
-- ref_idx = location in reference genome this is happening
-- mutation_type: 1 for delete, 2 for insert, 3 for SNP
-- insert_str:
-- insert_length:
-- snp_type: the nucleotid type this base changes into
CREATE TABLE IF NOT EXISTS aligned_bases (
    id SERIAL,
    ref_idx         INTEGER,
    mutation_type   INTEGER,
    insert_str      VARCHAR(50),
    insert_length   INTEGER,
    snp_type        INTEGER
);
