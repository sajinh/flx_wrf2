function mme_out_description(scheme_name,cabinet, hindcast_years, \
                           start_year, start_mon, season, lead_time, models)
begin
   
    mme_description                = scheme_name
    mme_description@source         = "APCC seasonal "+ cabinet +" based on "+mme_description+" method"
    hindcast_years@join            = ", "
    mme_description@hindcast_years = array_to_string(hindcast_years)
    models@join                    = " "
    mme_description@models_used     = array_to_string(models)
    forecast_info                  = lead_time +"-MON "+cabinet + " for "

    if str_upper(season) .eq. "NONE"
      forecast_info                = forecast_info + start_mon
    else
      forecast_info                = forecast_info + season
    end if
    
    if str_upper(cabinet) .eq. "HINDCAST"
      mme_description@hindcast_info  = forecast_info+" "+start_year
    else
      mme_description@forecast_info  = forecast_info+" "+start_year
    end if
  return(mme_description)
end
