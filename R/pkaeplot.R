#' Overlay adverse event duration on concentration-time curve plot
#'
#' \code{pkaeplot} creates a plot showing duration of adverse events in each subjects
#' overlayed on individual concentration-time curve
#'
#' The function returns a list of ggplot objects, including (1) plot for active arm,
#' (2) plot for placebo arm, and (3) legend.
#'
#'
#' @export
#' @param pk dataset containing individual pharmacokinetic profiles
#' @param ae dataset containing adverse event profiles
#' @param ae.data.first.day assign 1 if the first study day is recorded as
#' DAY=1 in AE dataset and DAY=0 in PK dataset
#' @param scale.y.log10 logical value controlling Y axis scale
#' @param x.range range of X axis to display
#' @param y.range range of Y axis to display for active arm
#' @param ae.col.var variable to set colors of AE (e.g. AE grade, AE type)
#' @param ae.col.name name appearing in legend for different colors of AE curves
#' @param ae.palette color palette for AE scales
#' @param pk.col.var variable to set colors of background PK curves (e.g. DOSE)
#' @param ggtheme specify ggplot theme (set NULL to use current theme)
#' @return list of ggplot objects
#' @examples
#' pkaeplot(pk = pk.sim.profile.sparse,
#'          ae = ae.sim.profile,
#'          ae.data.first.day = 1,
#'          y.range=c(0.3,3))
#'
#'
pkaeplot <-function(pk, ae, subj,
                    ae.data.first.day=0, scale.y.log10=T,
                    x.range=NULL, y.range=NULL,
                    ae.col.var="AE",
                    ae.col.name=NULL,
                    pk.col.var=NULL,
                    ae.palette=c("#56B4E9","#0072B2","#D55E00"),
                    ggtheme=NULL){

  #### Input data cleanup ####
  ae <- ae %>%
    mutate_(ae.col=ae.col.var) %>%
    mutate(ae.col=as.factor(ae.col)) %>%
    arrange(ae.col) %>%
    # Correct first day for AE dataset
    mutate(AESTDY = AESTDY-ae.data.first.day,
           AEEDDY = AEEDDY-ae.data.first.day) %>%
    select(ID,AE,ae.col,AESTDY,AEEDDY) %>%
    mutate(id.ae=1:nrow(.))

  ae.col.lev <- levels(ae$ae.col)

  pk <- pk %>%
    select(ID,DAY,DV) %>%
    mutate(no.pk=F)

  if(scale.y.log10) pk <- pk %>% filter(DV>0)

  if(is.null(ae.col.name)) ae.col.name <- ae.col.var
  if(is.null(ggtheme)) ggtheme <- theme_get()

  #### Prepare mock PK data for placebo group ####
  pk.nopk <-
    anti_join(subj, pk,by="ID") %>%
    mutate(min=min(pk$DAY),max=max(pk$DAY)) %>%
    arrange(ID) %>%
    # Give DV values in the increment of 1
    mutate(DV=1:nrow(.)) %>%
    gather(minmax,DAY,min,max) %>%
    select(ID,DAY,DV) %>%
    mutate(no.pk=T)

  n.nopk <- max(pk.nopk$DV)

  pk <- bind_rows(pk,pk.nopk)

  #### Combine PK and subject data ####
  pk <- full_join(pk,subj,by="ID")

  #### Assign colors for AE if specified ####
  if (is.null(pk.col.var)==F) {
    pk <-
      mutate_(pk,pk.col=pk.col.var) %>%
      mutate(pk.col=as.factor(pk.col))
  }

  #### Combine PK and AE data ####

  # Extrapolation of PK
  # Important if using observed PK profiles

  ## Select days that need extrapolation
  day.to.extrapolate <- ae %>%
    gather("STED", "DAY", AESTDY, AEEDDY) %>%
    select(ID,DAY) %>%
    # Choose AE start/end date not in PK dataset
    anti_join(pk,by=c("ID","DAY")) %>%
    # Remove subjects who had no PK record
    semi_join(pk,by="ID") %>%
    # Remove DAY==NA (e.g. unresolved AE)
    filter(is.na(DAY)==F)

  ## Remove days before/after PK data
  day.to.extrapolate <-
    pk %>% group_by(ID) %>%
    summarize(min=min(DAY),
              max=max(DAY)) %>%
    left_join(day.to.extrapolate,.,by="ID") %>%
    filter(DAY>=min&DAY<=max) %>%
    select(ID,DAY)

  ## Extrapolation calculation
  pk.to.combine <- pk

  if(nrow(day.to.extrapolate)>0){
    for(k in 1:nrow(day.to.extrapolate)){

      pk.each.id <- pk %>%
        filter(ID==day.to.extrapolate[[k,1]])

      if(scale.y.log10) pk.each.id$DV <- log(pk.each.id$DV)
      dv.ext <- approx(x=pk.each.id$DAY,
                       y=pk.each.id$DV,
                       xout=day.to.extrapolate[[k,2]],
                       rule=2)$y
      if(scale.y.log10) dv.ext <- exp(dv.ext)

      pk.to.combine <-
        bind_rows(pk.to.combine,
                  data.frame(day.to.extrapolate[k,],
                             DV=dv.ext,
                             no.pk=pk.each.id$no.pk[1])
        )
    }
  }


  # Subset of data for PK during AE
  pk.ae <-
    inner_join(pk.to.combine,ae,by="ID") %>%
    filter((DAY>=AESTDY|is.na(AESTDY)) & (DAY<=AEEDDY|is.na(AEEDDY)))

  # Subset of data for beginning and end of AE
  pk.ae.st.ed <-
    pk.ae %>%
    filter(DAY==AESTDY|DAY==AEEDDY)

  #### Set range of X and Y ####
  if (is.null(x.range)) x.range <- c(min(pk$DAY),max(pk$DAY))
  if (is.null(y.range)){
    pk2 <- filter(pk,no.pk==F)
    y.range <- c(min(pk2$DV),max(pk2$DV))
  }


  #### Function to plot PK and AE curve ####
  plot.pkae <-
    function(pk.plot,          # PK for all data points
             pk.ae.plot,       # PK for AE duration
             pk.ae.st.ed.plot, # PK for AE start and end day
             x.range,
             y.range,
             ae.palette, ae.col.lev, # Set colors for AE curves
             pk.col.var=NULL # Set colors for background PK curves
             ){


      # All PK curves
      if (is.null(pk.col.var)){
        g <-
          ggplot(pk.plot,aes(DAY, DV, group=id.ae)) +
          geom_line(aes(group=ID), color="grey50", alpha=0.4)
      } else {
        g <-
          ggplot(pk.plot,aes(DAY, DV, group=id.ae)) +
          geom_line(aes(group=ID, color=pk.col), alpha=0.4, show.legend = F)
      }

      # AE curves
      for (k in seq_along(ae.col.lev)){
        g <- g +
          geom_line(data=pk.ae.plot %>% filter(ae.col==ae.col.lev[k]),
                    linetype="solid", size=1.2, color = ae.palette[k]) +
          geom_point(data=pk.ae.st.ed.plot %>% filter(ae.col==ae.col.lev[k]),
                     size=2, color = ae.palette[k])
      }

      # Formatting
      g <- g +
        coord_cartesian(xlim = x.range, ylim = y.range) +
        ggtheme #+ theme(legend.position="none")

      return(g)
    }

  #### Function to set graph attributes ####

  set.plot.attr <- function(g){
    g <- g
    return(g)
  }

  #### Plot for subjects with PK ####
  g1 <-plot.pkae(pk.plot      = filter(pk,    no.pk==F),
                 pk.ae.plot   = filter(pk.ae, no.pk==F),
                 pk.ae.st.ed.plot = filter(pk.ae.st.ed,no.pk==F),
                 x.range,
                 y.range,
                 ae.palette, ae.col.lev,
                 pk.col.var)

  if(scale.y.log10) g1 <- g1+scale_y_log10()


  #### Plot for subjects without PK ####

  g2 <-plot.pkae(pk.plot      = filter(pk,    no.pk==T),
                 pk.ae.plot   = filter(pk.ae, no.pk==T),
                 pk.ae.st.ed.plot = filter(pk.ae.st.ed,no.pk==T),
                 x.range,
                 y.range = c(1,n.nopk),
                 ae.palette, ae.col.lev)


  #### Plot for legend ####

  g.for.legend <-
    ggplot(pk.ae,aes_string("DAY", "DV", group="id.ae", color="ae.col")) +
    geom_line(linetype="solid", size=1.2) +
    ggtheme + theme(legend.position="right") +
    scale_colour_manual(values=ae.palette,name=ae.col.name)

  # Extract legend
  tmp <- ggplot_gtable(ggplot_build(g.for.legend))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  g3 <- tmp$grobs[[leg]]


  #### Return value ####

  glist <- list(g1,g2,g3)

  return(glist)

}



#' Combine a list of figures from pkaeplot function
#'
#' \code{combine_pkaeplot} combines a list of figures from pkaeplot function
#'
#' Instead of using this function, you can write your own function to combine
#' list of figures for greater flexibility
#'
#' @export
#' @param glist list of figures generated with pkaeplot function
#' @param rel_heights numeric vector of relative columns heights of figures
#' for active and placebo arm (see cowplot::plot_grid for detail)
#' @param rel_widths numeric vector of relative columns widths of figures
#' and legend (see cowplot::plot_grid for detail)
#' @param legend whether to include legend in the plot
#' @return A graphic object containing figures for active and placebo arm
#' @examples
#' combine_pkaeplot(glist)
#'
#'
combine_pkaeplot <- function(glist,rel_heights=c(3,1),rel_widths=c(8,2),legend=TRUE){

  glist[[1]] <- glist[[1]] +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank())
  glist[[2]] <- glist[[2]] +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())

  plot <- cowplot::plot_grid(glist[[1]],glist[[2]],ncol=1,align="v",rel_heights = rel_heights)
  if(legend) plot <- cowplot::plot_grid(plot,glist[[3]],rel_widths = rel_widths)

  return(plot)
}
