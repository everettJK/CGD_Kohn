v <- toupper(readLines('vectorSeq'))
v <- gsub('\\s', '', v)
v <- gsub('\\d', '', v)
write(paste(v, collapse = ''), file = 'vectorSeq2')