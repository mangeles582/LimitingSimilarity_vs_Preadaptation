


##****************************************##
##* Maria Angeles Pérez-Navarro
##* King's College University
##* AlienImpacts project
##* plot required functions
##* Jan 2022
##* **************************************##

# ADam clarck to avoid zero floor in log transformation

log_zerosafe = function(x) {
  newx = x
  newx[x<0] = min(x[!is.na(x) & x>0])
  
  return(newx)
  
}


# variable_con <- "fun_w_dist"
# variable_cat <- "origin4"
# model <- mod7_plot
# transformation <- NA
# dataset <- plot_dis_df

eff_size_transf <- function(variable_con, # variable continuous x variable
                            variable_cat, #model cat variable, same as real data grouping var
                            model,
                            transformation,
                            dataset){
  
  ls <- tryCatch(list(x=seq(from=round(min(dataset[, variable_con], na.rm=T), 2),
                              to=round(max(dataset[, variable_con], na.rm=T),2), 
                              length.out = 100)), error = function(e) 60)
  names(ls) <- variable_con
  ef_var <- effects::Effect(c(variable_con, variable_cat),
                            model,ls)# only fixed effects
  df_eff_var <- data.frame(ef_var)
  
  if(missing(transformation)) {
    df_eff_var$variable_trans <- df_eff_var[, 1]
    names(df_eff_var)[2] <- "group"
    names(df_eff_var)[1] <- "variable"
    
  } else if(is.na(transformation)) {
      df_eff_var$variable_trans <- df_eff_var[, 1]
      names(df_eff_var)[2] <- "group"
      names(df_eff_var)[1] <- "variable"
      
    } else if (transformation=="scale"){
    
    df_eff_var$variable_trans <- scale(df_eff_var[, 1])
    names(df_eff_var)[2] <- "group"
    names(df_eff_var)[1] <- "variable"
    
  }else if(transformation=="log"){
    df_eff_var$variable_trans <- log(df_eff_var[, 1])
    names(df_eff_var)[2] <- "group"
    names(df_eff_var)[1] <- "variable"
    
  }else if(transformation=="log+1"){
    df_eff_var$variable_trans <- log(df_eff_var[, 1]+1)
    names(df_eff_var)[2] <- "group"
    names(df_eff_var)[1] <- "variable"
  }
  
  
  
  return(df_eff_var)
  
}



trend_plot <- function(real_data_x, real_data_y,#variable x and y real points include transformation if required
                       model, variable_con, # variable continuous x variable
                       variable_cat, #model cat variable, same as real data grouping var
                       transformation,
                       real_data_group,# real data grouping variable
                       mod_trend_group,# real data grouping variable
                       range_x0, range_x1, # xlim
                       range_y0, range_y1, # ylim
                       xtitle=NULL, # x variable name
                       ytitle="log(Cover percentage)", # y variable name
                       pvalue,# as value
                       model_r2,# as tibble
                       panel_title, # ex. a letter corresponding to panel order
                       color1="#EEC141",
                       color2="#494884",
                       text_x=1,
                       text_y=1,
                       legend_position="none"
                       ){
  
  ef_var <- effects::Effect(c(variable_con, variable_cat),
                            model,50)# only fixed effects
  df_eff_var <- data.frame(ef_var)
  
  if(missing(transformation)) {
    names(df_eff_var)[1] <- "variable"
    names(df_eff_var)[2] <- "group"
  } else if (transformation=="scale"){
    
    df_eff_var$variable <- scale(df_eff_var[, 1])
    names(df_eff_var)[2] <- "group"
    
  }else if(transformation=="log"){
    df_eff_var$variable <- log(df_eff_var[, 1])
    names(df_eff_var)[2] <- "group"
  }
  
  
  model_pv<- car::Anova(model)%>%
    as.data.frame()%>%
    rownames_to_column("variable")%>%
    rename(pv="Pr(>Chisq)")%>%
    filter(grepl(variable_con, variable))%>%
    filter(grepl(":",variable))%>%
    dplyr::select(pv)%>%
    pull()
  
  model_r2 <- MuMIn::r.squaredGLMM(model)%>%as.tibble()
  model_r2m <- model_r2%>%dplyr::select(R2m)%>%pull()
  model_r2c <- model_r2%>%dplyr::select(R2c)%>%pull()
  
  
  if(missing(range_x0)) {
    range_x0 <- min(c(real_data_x, df_eff_var$variable))
  } else {
    range_x0 <- range_x0
  }

  if(missing(range_x1)) {
    range_x1 <- max(c(real_data_x, df_eff_var$variable))
  } else {
    range_x1 <- range_x1
  }

  if(missing(range_y0)) {
    range_y0 <- min(real_data_y)
  } else {
    range_y0 <- range_y0
  }

  if(missing(range_y1)) {
    range_y1 <- max(real_data_y)
  } else {
    range_y1 <- range_y1
  }
  # 
  # 
  # if(missing(c(color1,color2))) {
  #   color1 <- "#EEC141"
  #   color2 <- "#494884"
  # } else {
  #   color1 <- color1
  #   color2 <- color2
  # }
  # 
  # if(missing(xtitle)) {
  #   xtitle <- NULL
  # } else {
  #   xtitle <- xtitle
  # }
  # 
  # if(missing(xtitle)) {
  #   ytitle <- NULL
  # } else {
  #   ytitle <- ytitle
  # }
  # 
  # 
  # if(missing(panel_title)) {
  #   panel_title <- NULL
  # } else {
  #   panel_title <- panel_title
  # }
  
  model_r2m <- model_r2$R2m
  model_r2c <- model_r2$R2c
  
  gg_trendplot <- ggplot()+
     geom_point(aes(x=real_data_x,
                    y=real_data_y,
                    color=real_data_group,
                    alpha=real_data_group),
                stroke=0.01)+
     scale_alpha_discrete(range = c(0.3, 0.05))+
     geom_line(aes(x=df_eff_var$variable,
                   y= df_eff_var$fit, 
                   color=df_eff_var$group),
               size=0.9,alpha=1)+
     geom_ribbon(aes(x=df_eff_var$variable,
                     ymin=df_eff_var$lower, 
                     ymax=df_eff_var$upper,
                     fill=df_eff_var$group),
                 alpha=0.3)+
     scale_color_manual(values=c(color1,
                                 color2))+
     scale_fill_manual(values=c(color1,
                                color2))+
     ylab ( ytitle) +
     xlab ( xtitle) +
     xlim(range_x0, range_x1)+
     ylim(range_y0, range_y1)+
     annotate("text", x=text_x, y=text_y, 
              label=paste0("R2m=", signif(model_r2m,3)), 
              hjust=0, size=2.6, family="serif")+
     annotate("text", x=text_x, y=text_y-0.3, 
              label=paste0("R2c=", signif(model_r2c,3)), 
              hjust=0, size=2.6, family="serif")+
     annotate("text", x=text_x, y=text_y-0.6, 
               label=paste0("Pv=", signif(model_pv,3)), 
               hjust=0, size=2.6, family="serif")+
    guides(color=F, size=F,alpha=F,
           fill = guide_legend(override.aes = list(alpha=1)))+
    labs(fill = "Origin")+
     ggtitle(panel_title)+
     theme_bw()+
     theme(legend.position=legend_position,
           legend.title= element_text(size=8),
           legend.text = element_text(size=7),
           axis.text.x = element_text(size=8, color = "black"),
           axis.text.y = element_text(size=8, color = "black"),
           axis.title.x = element_text(color="black", size=10 
           ),
           axis.title.y = element_text(
             color="black", size=10),
           text=element_text(family="serif"),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank())
  
  
  return(gg_trendplot)
  #print(gg_trendplot)
  
}



scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}





CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)




facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}



