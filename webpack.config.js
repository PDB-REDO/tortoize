const path = require('path');
const glob = require('glob-all');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const { PurgeCSSPlugin } = require("purgecss-webpack-plugin");
var TerserPlugin = require("terser-webpack-plugin");

const SCRIPTS = path.resolve(__dirname, "webapp");
const IMAGES = path.resolve(__dirname, "images");
const DEST = path.resolve(__dirname, "docroot");

module.exports = (env) => {

	const PRODUCTION = env != null && env.PRODUCTION;

	const webpackConf = {
		entry: {
			'index': path.resolve(SCRIPTS, 'index.js'),
			'PDB_logo_rect_medium': path.resolve(IMAGES, 'PDB_logo_rect_medium.svg')
		},

		output: {
			path: DEST,
			crossOriginLoading: 'anonymous',
			filename: "scripts/[name].js",
			clean: {
				// dry: true,
				keep(asset) {
					return asset.includes(".html")
				}
			}
		},

		module: {
			rules: [
				{
					test: /\.js$/,
					exclude: /node_modules/,
					use: {
						loader: "babel-loader",
						options: {
							presets: ['@babel/preset-env']
						}
					}
				},

				{
					test: /\.(sa|sc|c)ss$/i,
					use: [
						MiniCssExtractPlugin.loader,
						"css-loader",
						{
							loader: "sass-loader",
							options: {
								sassOptions: {
									silenceDeprecations: [
										"color-functions",
										"global-builtin",
										"import",
										"if-function"
									]
								}
							}
						}
					]
				},

				{
					test: /\.(woff(2)?|ttf)(\?v=[0-9]\.[0-9]\.[0-9])?$/,
					include: [
						path.resolve(__dirname, './node_modules/bootstrap-icons/font/fonts'),
						path.resolve(SCRIPTS, 'fonts')
					],
					type: 'asset/resource',
					generator: {
						filename: 'fonts/[name][ext]'
					}
				},

				{
					test: /\.svg$/,
					type: 'asset/resource',
					generator: {
						filename: 'images/[name][ext]'
					}
				}
			]
		},

		resolve: {
			extensions: ['.js', '.css', '.scss'],
		},

		optimization: {
			minimize: true,
			minimizer: [new TerserPlugin()]
		},

		target: 'web',

		plugins: [
			new MiniCssExtractPlugin({
				filename: "css/[name].css"
			}),
			new PurgeCSSPlugin({
				paths: glob.sync([
					`${SCRIPTS}/**/*`,
					`${DEST}/**/*`
				], { nodir: true })
			})
		]
	};

	if (PRODUCTION) {
		webpackConf.mode = "production";
	} else {
		webpackConf.mode = "development";
		webpackConf.devtool = 'source-map';
	}

	return webpackConf;
};
