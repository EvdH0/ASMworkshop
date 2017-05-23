#!/bin/bash

# Modified from https://stackoverflow.com/questions/28919191/write-script-to-create-multiple-users-with-pre-defined-passwords/28919425#28919425
# NOTE: Be sure to run this script with `sudo`.

# Read user and password
while read iuser ipasswd; do

  # Just print this for debugging.
  printf "\tCreating user: %s with password: %s\n" $iuser $ipasswd

  # Create the user with adduser (you can add whichever option you like).
  useradd -m -g workshop -s /bin/bash $iuser

  # Assign the password to the user.
  # Password is passed via stdin, *twice* (for confirmation).
  passwd $iuser <<< "$ipasswd"$'\n'"$ipasswd"

done < <(paste users.txt passwords.txt)
