function afs_getenv(env_name)
begin
  split_by_comma_kinds=(/"MME_MODELS","TRAINING_YEARS" /)
  if any(split_by_comma_kinds .eq. env_name)
    return(str_split(getenv(env_name),","))
  end if

  toint_kinds=(/"LEAD_TIME","START_YEAR"/)
  if any(toint_kinds.eq. env_name)
    return(stringtoint(getenv(env_name)))
  end if
   
  boolean_kinds=(/"AFS_Overwrite_Output"/)
  if any(boolean_kinds.eq. env_name)
    if ismissing(getenv(env_name))
      return(False)
    else
      value=(getenv(env_name))
      if str_upper(value).eq."FALSE"
        return(False)
      else
        return(True)
      end if
    end if
  end if
  return(getenv(env_name))
end

end
