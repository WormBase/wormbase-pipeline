environments{
   development{
   }
   test{
   }
   production{
      dataSource {
         dbCreate = "update" // one of 'create', 'create-drop', 'update', 'validate', ''
         // previous WBPS apollo used the database user 'apollo', in 2.6.1 it appears to be "apollo2_user"
         //username = "apollo"
         //password = "apollo"
         username = "apollo2_user"
         password = "**REDACTED**"
         driverClassName = "org.postgresql.Driver"
         dialect = "org.hibernate.dialect.PostgresPlusDialect"
         // previous WBPS apollo used the database name 'apollo', default in 2.6.1 is 'apollo-release-production'
         //url = "jdbc:postgresql://127.0.0.1/apollo"
         url = "jdbc:postgresql://127.0.0.1/apollo-release-production"
         properties {
            // See http://grails.org/doc/latest/guide/conf.html#dataSource for documentation
            jmxEnabled                    = false
            initialSize                   = 5
            maxActive                     = 50
            minIdle                       = 5
            maxIdle                       = 25
            maxWait                       = 10000
            maxAge                        = 10 * 60000
            timeBetweenEvictionRunsMillis = 5000
            minEvictableIdleTimeMillis    = 60000
            validationQuery               = "SELECT 1"
            validationQueryTimeout        = 3
            validationInterval            = 15000
            testOnBorrow = true
            testWhileIdle = true
            testOnReturn = false
            jdbcInterceptors = "ConnectionState"
            defaultTransactionIsolation = java.sql.Connection.TRANSACTION_READ_COMMITTED
         }
      }
      dataSource_chado {
         dbCreate = "update"
         username = "apollo2_user"
         password = "**REDACTED**"
      
         driverClassName = "org.postgresql.Driver"
         dialect = "org.hibernate.dialect.PostgresPlusDialect"
      
         // previous WBPS apollo used the database name 'chado', default in 2.6.1 is apollo-production-chado
         //url = "jdbc:postgresql://127.0.0.1/chado"
         url = "jdbc:postgresql://127.0.0.1/apollo-production-chado"
      
         properties {
            // See http://grails.org/doc/latest/guide/conf.html#dataSource for documentation
            jmxEnabled                    = false
            initialSize                   = 5
            maxActive                     = 50
            minIdle                       = 5
            maxIdle                       = 25
            maxWait                       = 10000
            maxAge                        = 10 * 60000
            timeBetweenEvictionRunsMillis = 5000
            minEvictableIdleTimeMillis    = 60000
            validationQuery               = "SELECT 1"
            validationQueryTimeout        = 3
            validationInterval            = 15000
            testOnBorrow = true
            testWhileIdle = true
            testOnReturn = false
            jdbcInterceptors = "ConnectionState"
            defaultTransactionIsolation = java.sql.Connection.TRANSACTION_READ_COMMITTED
         }
      }
   }
}


apollo{
   // sequence_search_tools section including in default 2.6.1 config, not sure if it should be used
   sequence_search_tools {
      blat_nuc {
         search_exe = "/data/blat/bin/blat"
         search_class = "org.bbop.apollo.sequence.search.blat.BlatCommandLineNucleotideToNucleotide"
         name = "Blat nucleotide"
         params = "-minScore=0 -minMatch=1"
         tmp_dir = "/data/blat/tmp"
      }
      blat_prot {
         search_exe = "/data/blat/bin/blat"
         search_class = "org.bbop.apollo.sequence.search.blat.BlatCommandLineProteinToNucleotide"
         name = "Blat protein"
         params = "-minScore=0 -minMatch=1"
         tmp_dir = "/data/blat/tmp"
      }
   }
   // rest of the apollo{} settings, below, copied from previous WBPS apollo
   default_minimum_intron_size = 1
   history_size = 0
   overlapper_class = "org.bbop.apollo.sequence.OrfOverlapper"
   use_cds_for_new_transcripts =  true
   feature_has_dbxrefs = false
   feature_has_attributes = false
   feature_has_pubmed_ids = false
   feature_has_go_ids = false
   feature_has_comments = true
   feature_has_status = true
   translation_table = "/config/translation_tables/ncbi_1_translation_table.txt"
   get_translation_code = 1
   common_data_directory = '/data/temp'

   // TODO: should come from config or via preferences database
   splice_donor_sites = ["GT"]
   splice_acceptor_sites = ["AG"]
   gff3.source = "iris"

   // TODO: uncomment this for production
//    google_analytics = System.getenv("WEBAPOLLO_GOOGLE_ANALYTICS_ID") ?: ["UA-62921593-1"]

   admin{
      username  = "**REDACTED**"
      password  = "**REDACTED**"
      firstName = "**REDACTED**"
      lastName  = "**REDACTED**"
   }
}


jbrowse {
    git {
        url = "https://github.com/GMOD/jbrowse"
        // updated to specify newer version of jbrowse
        // ATOW the latest release is 1.16.10-release, but at the time Apollo 2.6.1 was released the
        // latest was 1.16.9-release.   See https://github.com/GMOD/jbrowse/tags
        tag    = "1.16.9-release"
//      alwaysPull = true
//      alwaysRecheck = true
    }
   plugins {
      WebApollo{
         included = true
      }
         NeatHTMLFeatures{
         included = true
      }
         NeatCanvasFeatures{
         included = true
      }
         RegexSequenceSearch{
         included = true
      }
         HideTrackLabels{
         included = true
      }
   }
}
