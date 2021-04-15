getPacks <- function( packs ) {
  
  # Check wich packages are not intalled and install them
  if ( sum(!packs %in% installed.packages()[, 'Package'])) {
    install.packages( packs[ which(!packs %in% installed.packages()[, 'Package']) ], 
                      dependencies = T)
  }
  
  # Require all packages
  sapply(packs, require, character.only =  T)
  
}