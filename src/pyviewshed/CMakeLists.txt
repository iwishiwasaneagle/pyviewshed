set(LOCAL_SOURCES
        eventlist.cpp
        grass.cpp
        grid.cpp
        rbbst.cpp
        statusstructure.cpp
        utils.cpp
        viewshed.cpp
        visibility.cpp
        )

list(TRANSFORM LOCAL_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
target_sources(pyviewshed PRIVATE ${LOCAL_SOURCES})