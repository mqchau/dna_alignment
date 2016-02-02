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
-- mutation_type: 1 for delete, 2 for insert, 3 for SNP, 4 for match
-- insert_idx: start at 1, shows how far from the ref_idx this base is inserted in
-- new_base: of the SNP 
CREATE TABLE IF NOT EXISTS aligned_bases (
    id SERIAL,
    ref_idx         INTEGER,
    mutation_type   INTEGER,
    insert_idx      INTEGER,
    new_base        INTEGER
);


-- This table is used to store the result after we do local alignment after hashing match
-- ref_idx = same as aligned_bases
-- mutation_type: same as aligned_bases
-- insert_str: string that is inserted in here
-- new_base: of the SNP
CREATE TABLE IF NOT EXISTS mutation (
    id SERIAL,
    ref_idx         INTEGER,
    mutation_type   INTEGER,
    insert_str      VARCHAR(100),
    new_base        INTEGER
);
