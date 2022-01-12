
#pragma once

#include "pam_const.h"


class Notes {
public:

  struct Note {
    std::string key;
    std::string value;
  };

  std::vector<Note> notes;


  void add_note( std::string key , std::string value ) {
    notes.push_back( {key , value} );
  }

  std::string get_note( std::string key ) const {
    for (int i=0; i < notes.size(); i++) {
      if (key == notes[i].key) return notes[i].value;
    }
    return "";
  }

  bool note_exists( std::string key ) const {
    for (int i=0; i < notes.size(); i++) {
      if (key == notes[i].key) return true;
    }
    return false;
  }

  int get_num_notes() const {
    return notes.size();
  }

  void delete_note( std::string key ) {
    for (int i=0; i < notes.size(); i++) {
      if (key == notes[i].key) notes.erase(notes.begin()+i);
    }
  }
};


