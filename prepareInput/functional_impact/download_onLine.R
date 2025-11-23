rm(list = ls())
library(data.table)
library(dplyr)
source("functional_impact/functions.R")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                           DOWNLOAD
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
urls <- c('https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_8c4a9cd55012ec8625171b28bef8de51.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_6b12bdbb1049853d203dec2c91878642.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_ca679365690f60ec48cf6d3847a57111.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_96c5f5c8ead0f026ebc2e808e2881442.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_79a6964e7f70be49ae070f823e9fe390.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_c38258e13315e157b6de8200cb1bdd24.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_101975156b761098c50ab938a357e905.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_de065d3bbde1a339c71a0e9204a2b87d.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_de065d3bbde1a339c71a0e9204a2b87d.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_7f5f5b9a536aea4e3301425f3e2290a7.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_d4e181fe414cc5e4c1f7ace05ccd1754.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_ca61fb1df2e1a2ac1d288df4ae24f661.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_369e43b7fcd4e9a6aca9a3c61cbb8b0e.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_44ee6f30a6e9e2f8529a75f1827958fa.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_51b2f9d6dcc9dd4d9a2bb8f16513c8bd.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_96b17fca80863865b2a22ca3ae9e69b3.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_fb91f75d18573f2af1fef8e7673a757c.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_5995774e47cd1ad0879563b808492ca2.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_b5670a55f1ac9524bd6ec58e1200f39a.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_525d78f7891445cc8cc72c5d6ffab971.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_83576ad7fe7647aaf44872bafe37be7c.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_bb17dd770ae7c6c23833ec92b8bfb6f9.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_0308eb61a5509abff46690c87491dcca.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_a8f470cd86648df7d8443a04158bb39f.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_83a525023fea3ec6e7e93ba8cb5de853.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_4143364a882134eb05ad5ab9acd06e58.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_294467424ee49c7a80f252c60778f71e.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_d441280ce3ad5a219731387252a86d86.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_2aca601c5419920282f0b93f43ec27ba.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_5ed3b0fd156bc77956276486c511b1d2.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_13f7078eb256df9c61e0d54d0aae266b.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_97ae89a2d4e9ce830041e776f30eb8dc.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh38-v1.6_f20cdea6b2a97058f1a2548e68c8aae0.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_f10bdf8aa51d9d9082f16416a5f1d515.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_da86026d0b3a7b5d0ec300beee0ee067.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_89b04cf85b5e096e879be9eda0643418.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_5dbce014d66b706902c6228ebc44c055.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_1590dc46f4d62602833e4286dde255fa.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_0100bfc11cf43372e41e060dd003d798.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_d6c1aaebad6b64004aa299c2586b1860.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.6_1557e22a3279239947c898bae7dbddc7.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.7_da86026d0b3a7b5d0ec300beee0ee067.tsv.gz',
          'https://cadd.gs.washington.edu/static/finished/GRCh37-v1.7_f10bdf8aa51d9d9082f16416a5f1d515.tsv.gz')

path_save <- "../extdata/procInput/CADD_scores/onLine/"
dir.create(path_save, showWarnings = F, recursive = T)
file_names <- paste0(path_save, "file_", 1:length(urls), ".tsv.gz")
Map(function(u, f) download.file(u, f, mode = "wb"), urls, file_names)
