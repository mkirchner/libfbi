macro(MACRO_COPY_FILES GLOBPAT DESTINATION)
  file(GLOB COPY_FILES
    RELATIVE ${GLOBPAT}
    ${GLOBPAT}/*
    )
  add_custom_target(copy ALL
    COMMENT "Copying files: ${GLOBPAT}")

  foreach(FILENAME ${COPY_FILES})
    set(SRC "${GLOBPAT}/${FILENAME}")
    set(DST "${DESTINATION}/${FILENAME}"
    )

    add_custom_command(
      TARGET copy
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
      )
  endforeach(FILENAME)
endmacro(MACRO_COPY_FILES)
