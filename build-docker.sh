#!/bin/bash
#
# author: nils.hoffmann@isas.de
# date: August 2017
# 
NUMARGS=$#
DOCKER_NO_CACHE=""
for i in "$@"
do
case $i in
    -d=*|--directory=*)
    APP_DIR="${i#*=}"
    shift # past argument=value
    ;;
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
    --no-cache)
    DOCKER_NO_CACHE=" --no-cache"
    shift
    ;;
    *)
            # unknown option
    ;;
esac
done

CONFIG="./config.sh"

if [ -f "$CONFIG" ]; then
  source "$CONFIG"
else
  echo "File config.sh does not exist! Creating with default values!"
  echo -e "DOCKER_REGISTRY=\"do1-aps-feris.isas.de:5000\"" > "$CONFIG"
  echo -e "APP_IMAGE_BASE=\"isas\"" >> "$CONFIG"
  git add "$CONFIG"
fi

[[ -z "$DOCKER_REGISTRY" ]] && { echo -e "DOCKER_REGISTRY must be set in $CONFIG, example: DOCKER_REGISTRY=\"do1-aps-feris.isas.de:5000\""; exit 1; }
[[ -z "$APP_IMAGE_BASE" ]] && { echo -e "APP_IMAGE_BASE must be set in $CONFIG, example: APP_IMAGE_BASE=\"isas\""; exit 1; }

[[ -z "$APP_DIR" ]] && { echo "--directory= or -d= must be set!" ; exit 1; }
[[ -z "$APP_VERSION" ]] && { echo "--version= or -v= must be set!" ; exit 1; }

echo "Running $0"

function printHelp {
  echo -e "usage:"
  echo -e "\t-d=<APP_DIR> | --directory=<APP_DIR>>\t Use <APP_DIR> to find Dockerfile and application.sh in. (mandatory)"
  echo -e "\t-v=<VERSION> | --version=<VERSION>\t Use <VERSION> as the current version. (mandatory)"
  echo -e "\t--push \t Push the git tag and docker image to the remote repository and registry."
  echo -e "\t--tag  \t Tag the current build, both in git and in docker."
  echo -e "\t-? | --help \t Print this help."
  exit 1
}

[[ $NUMARGS -eq 0 || "$HELP" == true ]] && { printHelp; }


APP_CONFIG="$APP_DIR/application.sh"
APP_DOCKERFILE="$APP_DIR/Dockerfile.tmpl"

[[ ! -f "$APP_CONFIG" ]] && { echo "File $APP_CONFIG does not exist!"; exit 1; }
source "$APP_CONFIG"
[[ ! -f "$APP_DOCKERFILE" ]] && { echo "File $APP_DOCKERFILE does not exist!"; exit 1; }

echo "APP_DIR             = ${APP_DIR}"
echo "APP_VERSION         = ${APP_VERSION}"
echo "TAG (GIT + DOCKER)  = ${TAG}"
echo "PUSH (GIT + DOCKER) = ${PUSH}"
echo "APP_NAME            = ${APP_NAME}"
echo "APP_LICENSE         = ${APP_LICENSE}"
echo "APP_AUTHORS         = ${APP_AUTHORS}"
echo "APP_MAINTAINER      = ${APP_MAINTAINER}"
echo "APP_FROM_IMAGE      = ${APP_FROM_IMAGE}"
echo "APP_SCM             = ${APP_SCM}"


BUILD_TIMESTAMP="$(date --iso-8601='seconds')"
APP_LABEL_FILE="$APP_DIR/Dockerfile.labels"
# replace version
echo "Writing to $APP_LABEL_FILE"
echo "LABEL application.version=\"$APP_VERSION\"" > "$APP_LABEL_FILE"
echo "LABEL application.build.timestamp=\"$BUILD_TIMESTAMP\"" >> "$APP_LABEL_FILE"
echo "LABEL application.name=\"$APP_NAME\"" >> "$APP_LABEL_FILE"
echo "LABEL application.authors=\"$APP_AUTHORS\"" >> "$APP_LABEL_FILE"
echo "LABEL application.maintainer=\"$APP_MAINTAINER\"" >> "$APP_LABEL_FILE"
echo "LABEL application.license=\"$APP_LICENSE\"" >> "$APP_LABEL_FILE"
echo "LABEL application.from.image=\"$APP_FROM_IMAGE\"" >> "$APP_LABEL_FILE"
echo "LABEL application.scm=\"$APP_SCM\"" >> "$APP_LABEL_FILE"

echo "Contents of $APP_LABEL_FILE:"
cat "$APP_LABEL_FILE"
echo ""

APP_IMAGE_FILE="$APP_DIR/Dockerfile"

echo "Appending $APP_DOCKERFILE to $APP_IMAGE_FILE"
cat "$APP_DOCKERFILE" > "$APP_IMAGE_FILE"
echo "Replacing $APP_VERSION in $APP_IMAGE_FILE"
sed -ibak -e "s/\@APP_VERSION\@/${APP_VERSION}/g" "$APP_IMAGE_FILE"

echo "Appending $APP_LABEL_FILE to $APP_IMAGE_FILE"
cat "$APP_LABEL_FILE" >> "$APP_IMAGE_FILE"

echo "Contents of $APP_IMAGE_FILE:"

APP_IMAGE_TAG="$DOCKER_REGISTRY/$APP_IMAGE_BASE/$APP_NAME:$APP_VERSION"

echo "Adding and committing updated version file"
git add "$APP_LABEL_FILE" "$APP_IMAGE_FILE" "$APP_DOCKERFILE" && \
git commit -m "Updated version file of $APP_NAME to version $APP_VERSION."

if [ $TAG ]; 
then
    echo "Creating annotated git tag..."
    git tag -a "$APP_NAME-$APP_VERSION" -m "Release tag for $APP_NAME $APP_VERSION."
    [[ $? -ne 0 ]] && { echo "Git tag $APP_NAME-$APP_VERSION already exists! Please either remove it (git tag -d $APP_NAME-$APP_VERSION) or provide a higher version!" ; exit 1; }
fi

echo "Creating tagged docker images..."
echo "Building $APP_IMAGE_TAG..."
SRC_DIR=$PWD
cd "$APP_DIR"
docker build$DOCKER_NO_CACHE -t "$APP_IMAGE_TAG" .
[[ $? -ne 0 ]] && { echo "Docker build did not finish cleanly, please inspect output, correct and retry!"; exit 1;}
cd "$SRC_DIR"
if [ $PUSH ];
then
    echo "Pushing to git..."
    git push
    echo "Pushing to docker registry..."
    docker push "$APP_IMAGE_TAG"
fi
[[ $? -ne 0 ]] && { echo "Docker push of $APP_IMAGE_TAG did not finish cleanly, please inspect output, correct and retry!"; exit 1;}
echo "Done!"
