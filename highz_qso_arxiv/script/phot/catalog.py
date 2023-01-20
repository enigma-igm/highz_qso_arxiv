from wildhunt import catalog

if __name__ == '__main__':

    filename = '../../arxiv/hizqa_catalog.csv'
    ra_colname = 'ra'
    dec_colname = 'dec'

    # Instantiate a catalog from a csv file
    cat = catalog.Catalog('example', ra_column_name=ra_colname,
                              dec_column_name=dec_colname,
                              id_column_name='name',
                              datapath=filename)
                              
    # Do the online cross-match
    # cat.online_cross_match(survey='UKIDSSDR11LAS')
    # cat.online_cross_match(survey='VIKINGDR5')
    cat.online_cross_match(survey='DELS')