-- 1) Switch SQLite into tab-delimited output mode
.mode tabs

-- 2) Redirect all query output into a file named miRDB_subset.gmt
.once miRDB_subset.gmt

-- 3) Run the query:
SELECT 
    standard_name 'na',
    group_concat(symbol, '       ')
FROM gene_set gset
    INNER JOIN gene_set_gene_symbol gsgs ON gset.id = gsgs.gene_set_id
    INNER JOIN gene_symbol gsym ON gsym.id = gsgs.gene_symbol_id
WHERE collection_name = 'C3'
GROUP BY standard_name
ORDER BY standard_name ASC;

