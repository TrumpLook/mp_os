#include "cw/database/include/data_base.h"

int main()
{
    /*
     * конструктор с 1 параметром - режим in-memory
    */
    data_base db_1("NEW_DB", storage_interface<std::string, schemas_pool>::storage_strategy::filesystem);
    db_1.insert_schemas_pool("pl1", schemas_pool());
    db_1.insert_schema("pl1", "hr", schema());
    db_1.insert_table("pl1", "hr", "try", table());
    db_1.start_console_dialog();

    return 0;
}