# Git hooks for OpenGeoSys 6 #

This branch contains special git scripts which are executed on git events
like commit, push or pull.

Git looks for hooks in the .git/hooks directory within the work tree of a
local repository. Create a new local repository in this directory to manage
the hooks:

		$ cd .git/hooks
		$ git init
		$ cd ../..
		$ git fetch origin
		$ cd .git/hooks
		$ git pull .. remotes/origin/hooks
		$ cd ../..

