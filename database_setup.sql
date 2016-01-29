drop schema public cascade;
create schema public;

CREATE TABLE IF NOT EXISTS reference_hash (
    mer     VARCHAR(20),
    location INTEGER[],
    PRIMARY KEY (mer)
);
