
#pragma once

#include "pam_const.h"


class Notes {
public:

  struct Note {
    std::string key;
    std::string value;
  };

  std::vector<Note> notes;


  ~Notes() { finalize(); }


  void finalize() { notes = std::vector<Note>(); }


  void add_note( std::string key , std::string value ) {
    if ( key == "" ) return;
    int ind = get_note_index( key );
    if ( ind == -1 ) {
      notes.push_back( {key , value} );
    } else {
      notes[ind].value = value;
    }
  }


  void set_note( std::string key , std::string value ) {
    if ( key == "" ) return;
    add_note( key , value );
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


  int get_note_index( std::string key ) const {
    for (int i=0; i < notes.size(); i++) {
      if (key == notes[i].key) return i;
    }
    return -1;
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


