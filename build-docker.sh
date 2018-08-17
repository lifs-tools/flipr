#!/bin/bash
#
# author: nils.hoffmann@isas.de
# date: August 2017
# 
NUMARGS=$#
for i in "$@"
do
case $i in
    -v=*|--version=*)
    APP_VERSION="${i#*=}"
    shift # past argument=value
    ;;
    --push)
    PUSH=true
    shift # past argument with no value
    ;;
    --tag)
    TAG=true
    shift # past argument with no value
    ;;
    --help|-?)
    HELP=true
    shift # past argument with no value
    ;;
    *)
            # unknown option
    ;;
esac
done

source "./application.sh"
echo "Running $0"

function printHelp {
  echo -e "usage:"
  echo -e "\t-v=<VERSION> | --version=<VERSION>\t Use <VERSION> as the current version. (mandatory)"
  echo -e "\t--push \t Push the git tag and docker image to the remote repository and registry."
  echo -e "\t--tag  \t Tag the current build, both in git and in docker."
  echo -e "\t-? | --help \t Print this help."
  exit 1
}

[[ $NUMARGS -eq 0 || "$HELP" == true ]] && { printHelp; }

[[ -z "$APP_VERSION" ]] && { echo "--version= or -v= must be set!" ; exit 1; }

echo "APP_VERSION         = ${APP_VERSION}"
echo "TAG (GIT + DOCKER)  = ${TAG}"
echo "PUSH (GIT + DOCKER) = ${PUSH}"
echo "APP_NAME            = ${APP_NAME}"
echo "APP_LICENSE         = ${APP_LICENSE}"
echo "APP_AUTHORS         = ${APP_AUTHORS}"

BUILD_TIMESTAMP="$(date --iso-8601='seconds')"
APP_VERSION_FILE="version.R"
# replace version
echo "Writing to $APP_VERSION_FILE"
echo "application.version = \"$APP_VERSION\"" > "$APP_VERSION_FILE"
echo "application.build.timestamp = \"$BUILD_TIMESTAMP\"" >> "$APP_VERSION_FILE"
echo "application.name= \"$APP_NAME\"" >> "$APP_VERSION_FILE"
echo "application.authors = \"$APP_AUTHORS\"" >> "$APP_VERSION_FILE"
echo "application.license = \"$APP_LICENSE\"" >> "$APP_VERSION_FILE"

echo "Contents of $APP_VERSION_FILE:"
cat "$APP_VERSION_FILE"
echo ""

DOCKER_REGISTRY="do1-aps-feris.isas.de:5000"
APP_IMAGE_BASE="isas"

# suprg code + shiny application
APP_IMAGE_TAG="$DOCKER_REGISTRY/$APP_IMAGE_BASE/$APP_NAME:$APP_VERSION"
APP_IMAGE_FILE="docker/suprg/Dockerfile"

echo "Adding and committing updated version file"
git add "$APP_VERSION_FILE" && \
git commit -m "Updated version file of $APP_NAME to version $APP_VERSION."

if [ $TAG ]; 
then
    echo "Creating annotated git tag..."
    git tag -a "$APP_NAME-$APP_VERSION" -m "Release tag for $APP_NAME $APP_VERSION."
    [[ $? -ne 0 ]] && { echo "Git tag $APP_NAME-$APP_VERSION already exists! Please either remove it or provide a higher version!" ; exit 1; }
fi

echo "Creating tagged docker images..."
echo "Building $APP_IMAGE_TAG..."
docker build --file "$APP_IMAGE_FILE" -t "$APP_IMAGE_TAG" .

if [ $PUSH ];
then
    echo "Pushing to git..."
    git push
    echo "Pushing to docker registry..."
    docker push "$APP_IMAGE_TAG"
fi
echo "Done!"
