odbc <- function(dsn, sqtable, fields)
{
  package.already.there <- any(search()=="package:RODBC")
  if(!package.already.there)
    require(RODBC, quietly=TRUE)
  #specify full filepath for database!
  channel <- odbcConnectAccess2007(database)
  query <- paste("SELECT ", fields, " FROM [", sqtable, "]", sep="")
  imported <- sqlQuery(channel, query)
  odbcClose(channel)
  if(!package.already.there)
    detach("package:RODBC")
  return(imported)
}
